#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis code for article: Image-Based & Machine Learning-Guided Multiplexed Serology Test for SARS-CoV-2

Methods to generate per-well predictions with a given classifier

@author: Christian Guckelsberger, Lassi Paavolainen
"""
import re
import os
import pandas as pd
import numpy as np
import helpers
import gc
from pathlib import Path
from sklearn import preprocessing

#
# Calculate predictions individually for each antibody and well
#
def predict_per_well(classifier, classes, df_features_meta, dilution, channel):
    # Picking well identifiers from df_features_meta: only those wells that have features assigned
    df_agg = pd.DataFrame({'well': df_features_meta.well.unique()})

    print("\nCalculating predictions for antibody {} and dilution {}".format(channel, dilution))

    # Subset features for prediction on this channel 
    features = [col for col in df_features_meta.columns if (("NUCLEUS" in col or "CELL" in col or "DONUT" in col) and not 'Alexa' in col and not 'NUCLEUS.NUMBER.OF.COMPONENTS' in col and not ('METCENTER' in col))]
    features = features + [col for col in df_features_meta.columns if channel in col]
    features.sort()
    
    col_name_cell_pred = 'Prediction {}'.format(channel)
    col_names = dict()
    for i, c in classes.items():
        col_names[i] = 'Ratio {} {}'.format(c, channel)

    # Perform prediction for antibody-specific data
    print("Calculating per-cell predictions")
    df_features_meta[col_name_cell_pred] = classifier.predict(df_features_meta[features])
    
    # Aggregate per-well information
    print("Calculating per-well ratios")
    for ix, well in df_agg.well.items():
        features_per_well = df_features_meta[df_features_meta.well==well]
      
        n = features_per_well[col_name_cell_pred].count()
        df_agg.loc[ix, "Cell Count"] = n
        for i, c in classes.items():
            df_agg.loc[ix, col_names[i]] = (features_per_well[col_name_cell_pred]==i).sum() / n
            
    return df_agg, df_features_meta


def predict_plate(classifier, classifier_type, path_plate, path_output, antibody, classes, df_controls, df_meta, normalisation_type, dilution_antibody):

    # Plate id from filename
    pat_plate = re.compile(r'([SNRM]_\d+_\d+)')
    plate_name = path_plate.name.split('.')[0]
    plate_id = pat_plate.search(plate_name).groups()[0]    
    
    # Features to be discarded
    features_discard = ['NUCLEUS.ELLIPSE-ORIENT']

    # Quality control
    # Feature used for QC to check if controls are transferred correctly
    # And Std multiplier for QC feature. The range is mean +- qc_mult * std
    qc_feat = 'NUCLEUS.INTENSITY-MEAN.DAPI'
    qc_mult = 1.0

    # Only retain metadata for requested plate
    df_meta_plate = df_meta[df_meta['Destination Plate Barcode'] == plate_id]

    #
    # Load features for specific plate
    print('Loading feature data for plate {}'.format(path_plate))
    df_features_plate = pd.read_hdf(path_plate)

    # Drop irrelevant features (ignore: only drop available columns; don't break if one feature to discard is not in the data)
    df_features_plate.drop(features_discard, axis=1, inplace=True, errors='ignore')

    #
    # Construct overall data set
    df_features_meta = pd.merge(df_meta_plate, df_features_plate, how='inner', on='well')
    df_features_meta.dropna(inplace=True)

    # Add index column to uniquely identify individual cells later once the dataframe has been filtered
    df_features_meta['index'] = np.arange(0,df_features_meta.shape[0])

    #
    # Construct mask for features that are relevant to training.
    # Sort to make sure that the order is identical to the order of features in the trained model
    features = [c for c in df_features_plate.columns if (("NUCLEUS" in c or "CELL" in c or "DONUT" in c) and not ('NUCLEUS NUMBER OF COMPONENTS' in c) and not ('METCENTER' in c))]
    features.sort()

    del(df_features_plate)
    gc.collect()
    
    #
    # Filter data for antibody-specific dilution and normalise
    d = dilution_antibody[antibody]
    print("\nNormalising data for dilution {}".format(d))
    df_features_meta_dil = df_features_meta[df_features_meta.Dilution==d].copy()

    # Check if sufficient number of controls have been transferred
    controls_pos = df_controls[(df_controls.Dilution == d) & (df_controls.Type == 'positive')].Content.to_list()
    controls_pos_transferred = df_meta_plate.loc[df_meta_plate.Content.isin(controls_pos),['Transferred','well']]
    controls_neg = df_controls[(df_controls.Dilution == d) & (df_controls.Type == 'negative')].Content.to_list()
    controls_neg_transferred = df_meta_plate.loc[df_meta_plate.Content.isin(controls_neg),['Transferred','well']]
            
    # Fit normalisers on control wells for dilution (i.e. labelled and unlabelled samples)
    normalisation_controls_wells = helpers.select_normalisation_wells(normalisation_type, controls_pos_transferred, controls_neg_transferred)
    scaler = preprocessing.StandardScaler().fit(df_features_meta_dil[df_features_meta_dil.well.isin(normalisation_controls_wells)][features])

    df_features_meta_dil[features] = scaler.transform(df_features_meta_dil[features])
    df_features_meta_dil.drop(['CLASS.IgM','CLASS.IgG','CLASS.IgA'], axis=1, inplace=True, errors='ignore')
   
    #
    # Keep relevant, raw (unnormalised) intensity data to add as quality control to ratios later
    df_intensities = df_features_meta[[
        'index','well',
        'CELL.INTENSITY-MEAN.Alexa.IgM','CELL.INTENSITY-MEAN.Alexa.IgG','CELL.INTENSITY-MEAN.Alexa.IgA',
        'DONUT.INTENSITY-MEAN.Alexa.IgM','DONUT.INTENSITY-MEAN.Alexa.IgG','DONUT.INTENSITY-MEAN.Alexa.IgA',
        'NUCLEUS.INTENSITY-MEAN.DAPI'
    ]].copy()
    
    #
    # Calculate per-cell and per-well predictions for different dilutions
    df_pred, _ = predict_per_well(classifier, classes, df_features_meta_dil, d, antibody)

    # Merge in metadata, and remove columns that are not required for further use
    df_agg = df_pred.merge(df_meta_plate, how='left', on='well')
    df_agg.drop(['Sample_No', 'Layout name', 'Secondary'], axis=1, inplace=True, errors='ignore')

    #
    # Export to .csv
    Path(path_output).mkdir(parents=True, exist_ok=True)
    path_well_results = Path(path_output) / "ratios_per_well_plate_{}_ab_{}.csv".format(plate_id, antibody)
    print("\nExporting per-well prediction results to {}\n".format(path_well_results))
    df_agg.to_csv(path_well_results, index=False)  


def predict(classifier, classifier_type, path_data, path_output, antibody, classes, df_controls, df_meta, normalisation_type, dilution_antibody):

    path_list = Path(path_data).glob('**/*.hdf5')
    for filepath in path_list:
        predict_plate(classifier, classifier_type, filepath, path_output, antibody, classes, df_controls, df_meta, normalisation_type, dilution_antibody)
