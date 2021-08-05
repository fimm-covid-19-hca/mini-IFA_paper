#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Constructs per-cell training and training/validation data sets for regular training and cross-validation

@author: Christian Guckelsberger, Lassi Paavolainen
"""
import re
import helpers
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn import preprocessing
        
'''
Constructs data for leave-one-out cross-validation. 
Takes input feature data from path_data, and outputs the cross-validation files to path_output. These paths should be different in order to not mix up the two types of data.
'''
def construct_data_cross_validation(path_data, path_output, classes, df_controls, df_meta, normalisation_type, antibody, dilution_antibody):

    # Accumulate feature data from multiple plates and normalise using per-plate controls
    path_list = Path(path_data).glob('**/*.hdf5')
    gt_plate = dict()
    plate_ids = []
    for filepath in path_list:
        
        # Read plate 
        pat_plate = re.compile(r'([SNRM]_\d+_\d+)')
        plate_name = filepath.name.split('.')[0]
        plate_id = pat_plate.search(plate_name).groups()[0]
        plate_ids.append(plate_id)

        print("Reading plate {}".format(plate_id))
        df_features = pd.read_hdf(filepath)

        # Drop features for consistency
        df_features.drop(['NUCLEUS.ELLIPSE-ORIENT'],axis=1,inplace=True,errors='ignore')

        # Select metadata for this plate and merge into feature data
        df_meta_plate = df_meta[df_meta["Destination Plate Barcode"]==plate_name]
        df_features = pd.merge(df_features, df_meta_plate, how='left', on='well')

        # Subselect training features. Should have a label assigned to at least one channel. 
        df_features_labelled = df_features[(df_features['CLASS.IgM']+df_features['CLASS.IgG']+df_features['CLASS.IgA'])>0].copy()
        
        # Filter out rows with NaN's if any
        df_features_labelled.dropna(inplace=True)
        
        def sample_count(df, antibody):
            col_class = 'CLASS.{}'.format(antibody)
            labels, counts = np.unique(df[df[col_class]>0][col_class], return_counts=True)
            for i in range(0, len(counts)):
                print("%s: %d" % (classes[labels[i]], counts[i]))
            
        # Filtering for specific dilutions
        # Print info on class counts for each channel for all and specific dilution
        df_features_labelled_dil = df_features_labelled.copy()
        dil_antibody = dilution_antibody[antibody]
        col_class = 'CLASS.{}'.format(antibody)
        df_features_labelled_dil.loc[df_features_labelled_dil.Dilution!=dil_antibody, col_class] = 0

        df_features_labelled_dil = df_features_labelled_dil[(df_features_labelled_dil['CLASS.IgM']+df_features_labelled_dil['CLASS.IgG']+df_features_labelled_dil['CLASS.IgA'])>0]

        print("\nLabelled samples for antibody {} for dilution {}:".format(antibody, dil_antibody))
        sample_count(df_features_labelled_dil, antibody)

        features = [c for c in df_features.columns if ("NUCLEUS" in c or "CELL" in c or "DONUT" in c)]

        #
        # Normalise data for dilutions 25 and 100 separately  
        print("\nNormalising for dilution {}".format(dil_antibody))

        # Flag wells that do not fulfill QC check
        flagged_wells = helpers.well_qc(df_features, dilution=dil_antibody)

        # Check if sufficient number of controls have been transferred
        controls_pos = df_controls[(df_controls.Dilution == dil_antibody) & (df_controls.Type == 'positive')].Content.to_list()
        controls_pos_transferred = df_meta_plate.loc[df_meta_plate.Content.isin(controls_pos),['Transferred','well']]
        controls_neg = df_controls[(df_controls.Dilution == dil_antibody) & (df_controls.Type == 'negative')].Content.to_list()
        controls_neg_transferred = df_meta_plate.loc[df_meta_plate.Content.isin(controls_neg),['Transferred','well']]
        helpers.check_sufficient_transfer(controls_pos_transferred, controls_neg_transferred, flagged_wells)
            
        # Fit normalisers on control wells for dilution (i.e. labelled and unlabelled samples)
        normalisation_controls_wells = helpers.select_normalisation_wells(normalisation_type, controls_pos_transferred, controls_neg_transferred)
        df_features_normalisation = df_features[df_features['well'].isin(normalisation_controls_wells)]
        scaler = preprocessing.StandardScaler().fit(df_features_normalisation[features])
            
        # Select labelled data for specific dilution and normalise
        df_features_labelled_norm = df_features_labelled_dil[df_features_labelled_dil.Dilution==dil_antibody].copy()
        df_features_labelled_norm[features] = scaler.transform(df_features_labelled_norm[features])
            
        # Store for this specific plate and dilution
        gt_plate[plate_id] = df_features_labelled_norm
        
        print()

    # Construct and export leave-one-plate-out-cross-validation data
    # For N input data plates, we have N folds of N-1 training and 1 validation set
            
    print("Exporting cv data for dilution {} and antibody {}".format(str(dil_antibody), antibody)) 

    # Determine relevant featurees for channel and dilution: DAPI, per-channel ALEXA, class label
    # Also drop NUCLEUS.NUMBER.OF.COMPONENTS column as 0 for all samples
    features_channel = [col for col in df_features.columns if (('NUCLEUS' in col or 'CELL' in col or 'DONUT' in col) and not ('Alexa' in col) and not ('NUCLEUS.NUMBER.OF.COMPONENTS' in col))]
    features_channel = features_channel + [col for col in df_features.columns if str(antibody) in col]

    # Sort by columns alphabetically. To make sure that columns have same order everywhere (training, prediction)
    features_channel.sort()

    # Leave one plate out for validation, use others for training
    col_class = 'CLASS.{}'.format(antibody)
    for validation_plate in plate_ids: 
        
        df_validation_dil = gt_plate[validation_plate]
        df_validation_dil = df_validation_dil[df_validation_dil[col_class]!=0]
        
        df_train_dil = pd.DataFrame()
        train_plates = []
        for plate_id in gt_plate.keys():
            if plate_id != validation_plate:
                train_plates.append(plate_id)
                
                df_train_plate = gt_plate[plate_id]
                df_train_plate = df_train_plate[df_train_plate[col_class]!=0]
                df_train_dil = df_train_dil.append(df_train_plate)
        
        # Filter for relevant features for this channel
        df_validation_dil = df_validation_dil[features_channel]
        df_train_dil = df_train_dil[features_channel]

        # Split up data into classes and drop class label colum
        df_validation_dil_class = dict()
        df_train_dil_class = dict()
        for label in classes.keys():
            df_validation_dil_class[label] = df_validation_dil[df_validation_dil[col_class]==label].drop([col_class], axis=1)
            df_train_dil_class[label] = df_train_dil[df_train_dil[col_class]==label].drop([col_class], axis=1)

        # Combine classes into datasets
        X_validate, X_train = pd.DataFrame(), pd.DataFrame()
        y_validate, y_train = pd.DataFrame(), pd.DataFrame()
        for label in classes.keys():
            X_validate = pd.concat([X_validate, df_validation_dil_class[label]], axis=0, ignore_index=True)
            X_train = pd.concat([X_train, df_train_dil_class[label]], axis=0, ignore_index=True)
            y_validate = pd.concat([y_validate, pd.DataFrame({'Label':np.full(df_validation_dil_class[label].shape[0], label, dtype=np.int8)})])
            y_train = pd.concat([y_train, pd.DataFrame({'Label':np.full(df_train_dil_class[label].shape[0], label, dtype=np.int8)})])

        print("Leaving out plate {}. Size of train and validation sets: {}, {}".format(validation_plate, X_train.shape[0], X_validate.shape[0]))

        # Export data. Produce two hdf files (X,y) with 2 groups each (train, validate)
        Path(path_output).mkdir(parents=True, exist_ok=True)
        filename = path_output / 'leave_{}_out_train_on_{}_ab_{}_and_DAPI_dil_{}.hdf5'.format(validation_plate, "_and_".join(train_plates), antibody, dil_antibody)
        X_train.to_hdf(filename, key='X_train', mode='w')
        X_validate.to_hdf(filename, key='X_validate')
        y_train.to_hdf(filename, key='y_train')
        y_validate.to_hdf(filename, key='y_validate')


'''
Constructs a single training data set (no train/validation split)
Takes input feature data from path_data, and outputs the training data files to path_output. These paths should be different in order to not mix up the two types of data.
'''
def construct_training_data(path_data, path_output, classes, df_controls, df_meta, normalisation_type, antibody, dilution_antibody):

    # Accumulate feature data from multiple plates and normalise using per-plate controls
    path_list = Path(path_data).glob('**/*.hdf5')
    gt = pd.DataFrame()
    for filepath in path_list:
        
        # Read plate 
        pat_plate = re.compile(r'([SNRM]_\d+_\d+)')
        plate_name = filepath.name.split('.')[0]
        plate_id = pat_plate.search(plate_name).groups()[0]        
        
        print("Reading plate {}".format(plate_id))
        df_features = pd.read_hdf(filepath)

        # Drop features for consistency
        df_features.drop(['NUCLEUS.ELLIPSE-ORIENT'],axis=1,inplace=True,errors='ignore')

        # Select metadata for this plate and merge into feature data
        df_meta_plate = df_meta[df_meta["Destination Plate Barcode"]==plate_id]
        df_features = pd.merge(df_features, df_meta_plate, how='left', on='well')

        # Subselect training features. Should have a label assigned to at least one channel. 
        df_features_labelled = df_features[(df_features['CLASS.IgM']+df_features['CLASS.IgG']+df_features['CLASS.IgA'])>0].copy()

        # Filter out rows with NaN's if any
        df_features_labelled.dropna(inplace=True)

        def sample_count(df, antibody):
            col_class = 'CLASS.{}'.format(antibody)
            labels, counts = np.unique(df[df[col_class]>0][col_class], return_counts=True)
            for i in range(0, len(counts)):
                print("%s: %d" % (classes[labels[i]], counts[i]))
            
        # Filtering for specific dilutions
        # Print info on class counts for each channel for all and specific dilution
        df_features_labelled_dil = df_features_labelled.copy()
        dil_antibody = dilution_antibody[antibody]
        col_class = 'CLASS.{}'.format(antibody)
        df_features_labelled_dil.loc[df_features_labelled_dil.Dilution!=dil_antibody, col_class] = 0

        df_features_labelled_dil = df_features_labelled_dil[(df_features_labelled_dil['CLASS.IgM']+df_features_labelled_dil['CLASS.IgG']+df_features_labelled_dil['CLASS.IgA'])>0]

        print("\nLabelled samples for antibody {} for dilution {}:".format(antibody, dil_antibody))
        sample_count(df_features_labelled_dil, antibody)
        
        features = [c for c in df_features.columns if ("NUCLEUS" in c or "CELL" in c or "DONUT" in c)]

        # Normalise data for dilutions 25 and 100 separately  
        d = dilution_antibody[antibody]
        print("\nNormalising for dilution {}".format(d))
            
        # Flag wells that do not fulfill QC check
        flagged_wells = helpers.well_qc(df_features, dilution=d)

        # Check if sufficient number of controls have been transferred
        controls_pos = df_controls[(df_controls.Dilution == d) & (df_controls.Type == 'positive')].Content.to_list()
        controls_pos_transferred = df_meta_plate.loc[df_meta_plate.Content.isin(controls_pos),['Transferred','well']]
        controls_neg = df_controls[(df_controls.Dilution == d) & (df_controls.Type == 'negative')].Content.to_list()
        controls_neg_transferred = df_meta_plate.loc[df_meta_plate.Content.isin(controls_neg),['Transferred','well']]
        helpers.check_sufficient_transfer(controls_pos_transferred, controls_neg_transferred, flagged_wells)
                
        # Fit normalisers on control wells for dilution (i.e. labelled and unlabelled samples)
        normalisation_controls_wells = helpers.select_normalisation_wells(normalisation_type, controls_pos_transferred, controls_neg_transferred)
        df_features_normalisation = df_features[df_features['well'].isin(normalisation_controls_wells)]
        scaler = preprocessing.StandardScaler().fit(df_features_normalisation[features])
            
        # Select labelled data for specific dilution and normalise
        df_features_labelled_norm = df_features_labelled_dil[df_features_labelled_dil.Dilution==d].copy()
        df_features_labelled_norm[features] = scaler.transform(df_features_labelled_norm[features])

        # Append to existing data
        gt = gt.append(df_features_labelled_norm) 
        
        print()
            
    # Construct and export training data
    print("Exporting training data for dilution {} and antibody {}".format(str(dil_antibody), antibody)) 
    
    # Extract labelled data for channel and dilution: DAPI, per-channel ALEXA, class label
    # Also drop NUCLEUS.NUMBER.OF.COMPONENTS column as 0 for all samples
    features_channel = [col for col in df_features.columns if (('NUCLEUS' in col or 'CELL' in col or 'DONUT' in col) and not ('Alexa' in col) and not ('NUCLEUS.NUMBER.OF.COMPONENTS' in col))]
    features_channel = features_channel + [col for col in df_features.columns if str(antibody) in col]
    df_features_labels = gt[features_channel]

    # Split up data into classes and drop class label column
    col_class = 'CLASS.{}'.format(antibody)
    df_class = dict()
    for label in classes.keys():
        df_class[label] = df_features_labels[df_features_labels[col_class]==label].drop([col_class], axis=1)
    
    # Combine classes into datasets
    X = pd.DataFrame()
    y = pd.DataFrame()
    for label in classes.keys():
        X = pd.concat([X, df_class[label]], axis=0, ignore_index=True)
        y = pd.concat([y, pd.DataFrame({'Label':np.full(df_class[label].shape[0], label, dtype=np.int8)})])
        
    # Sort by columns alphabetically. To make sure that columns have same order everywhere (training, prediction)
    X.sort_index(axis=1, inplace=True)
    
    # Save to harddrive
    Path(path_output).mkdir(parents=True, exist_ok=True)
    filename = path_output / 'train_ab_{}_and_DAPI_dil_{}.hdf5'.format(antibody, dil_antibody)
    X.to_hdf(filename, key='X', mode='w')
    y.to_hdf(filename, key='y')

    print()
