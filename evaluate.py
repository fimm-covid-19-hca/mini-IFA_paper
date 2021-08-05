#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis code for article: Image-Based & Machine Learning-Guided Multiplexed Serology Test for SARS-CoV-2

Evaluation methods to compare per-well positivity predictions against ELISA and visual ground truth data

@author: Christian Guckelsberger, Lassi Paavolainen
"""
import os
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import confusion_matrix, make_scorer, recall_score
from sklearn.metrics import auc as calc_auc

'''
Metrics
'''
# The proportion of actual positives that are correctly identified as such.
# The extent to which actual positives are not overlooked (so false negatives are few).
def score_sensitivity(y_true, y_pred, C=None):
    if not (type(C) is np.ndarray and C.shape==(2,2)):
        C = confusion_matrix(y_true, y_pred)
    _, _, fn, tp = C.ravel()
    return 0 if tp==0 else tp / (tp+fn)

# The proportion of actual negatives that are correctly identified as such.
# The extent to which actual negatives are classified as such (so false positives are few).
def score_specificity(y_true, y_pred, C=None):
    if not (type(C) is np.ndarray and C.shape==(2,2)):
        C = confusion_matrix(y_true, y_pred)
    tn, fp, _, _ = C.ravel()
    return tn / (tn+fp)

# Composite measure of specificity and sensitivity which retains the original range between [0,1].
# SpecSens = Specificity x Sensitivity
def score_specsens(y_true, y_pred, C=None):
    if not (type(C) is np.ndarray and C.shape==(2,2)):
        C = confusion_matrix(y_true, y_pred)
    return (score_specificity(y_true, y_pred, C) * score_sensitivity(y_true, y_pred, C))


def setup_scores():
    # Define scores to report on
    scores = {'specificity':make_scorer(score_specificity),
            'sensitivity':make_scorer(score_sensitivity),
            'specsens':make_scorer(score_specsens)}

    # Define score for which the model should be tuned
    score_tuning = 'specsens'

    return (scores, score_tuning)


'''
Function to find thresholds for given metrics
Thus can generate data for ROC curve and AUC metric
'''
def threshold(y_gt, y_pred_ratio, scores):
    
    # Configure thresholds for ROC-curve 
    # Rather than using linear thresholding, increase threshold width exponentially
    # exp_multiplier determines strength of exponential growth
    exp_multiplier = 4
    thresholds_lin = np.linspace(0, 1, num=5001)
    thresholds = (np.exp(exp_multiplier*thresholds_lin)-1) * (1.0/(np.exp(exp_multiplier)-1))
    threshold_steps = thresholds.shape[0]

    scores_val = dict()
    for score in scores:
        scores_val[score] = np.zeros(threshold_steps)
    
    for i in range(0,threshold_steps):
        y_pred = (y_pred_ratio >= thresholds[i]).astype('int')
        
        # Only for binary classification. Important to specify labels in case ground truth and predictions are complete 0 or 1 vectors
        C = confusion_matrix(y_gt, y_pred, labels=[0,1])
        for score in scores:
            scores_val[score][i] = scores[score]._score_func(y_gt, y_pred, C)
    
    auc = calc_auc(scores_val['specificity'], scores_val['sensitivity'])
    df_scores = pd.DataFrame.from_dict(scores_val)
    df_scores["Threshold"] = thresholds 

    return auc, df_scores


def evaluate_visual_gt(path_data, path_visual_ground_truth, protein, antibody, dilution_antibody):
    
    col_label_predictions = "Ratio positive"

    scores, _ = setup_scores()

    #
    # Load visual ground truth data
    print('Loading visual ground truth data from {}'.format(path_visual_ground_truth))
    df_gt = pd.read_csv(path_visual_ground_truth)
    col_gt_binary = "positive binary wo unclear"
    df_gt = df_gt[['Sample Name','Sample',"antibody", "Destination Plate Barcode", col_gt_binary]]
    df_pred_agg = pd.DataFrame()
    
    print('Loading per-well predictions and comparing against visual ground truth')
    file_list = [fn for fn in os.listdir(Path(path_data)) if (fn.endswith("csv") and (antibody in fn))]
    for filename in file_list:
        file_path = Path(path_data) / filename

        # Load per-well predictions
        df_predictions = pd.read_csv(file_path)
                
        # Exclude wells that have not transferred
        wells_excluded = (df_predictions.Transferred==False).sum()
        if wells_excluded>0:
            print("Excluding {} wells that have not transferred from evaluation".format(wells_excluded))
        df_predictions = df_predictions[df_predictions.Transferred==True]
        
        pred_cols = ["well", "Destination Plate Barcode", "Sample Name", "Sample", "Dilution", "Transferred"] + [c for c in df_predictions.columns if ('positive' in c)]
        df_predictions = df_predictions[pred_cols]

        #
        # Calculate dataset for per-channel continuous predictions and binary ground truth        
        col_pred = "{} {}".format(col_label_predictions, antibody)

        # Only retain per-well ratios for wells that contain a sample for which ground truth exists explicitly
        df_pred_ch = df_predictions[df_predictions.Dilution == dilution_antibody[antibody]]
        df_gt_ch = df_gt[df_gt["antibody"]==antibody].drop("antibody", axis=1)
        df_pred_gt = df_pred_ch.merge(df_gt_ch, on=['Sample Name', "Sample", "Destination Plate Barcode"], how='inner')
        df_pred_gt = df_pred_gt[['Sample Name', 'Sample', 'Destination Plate Barcode', col_pred, col_gt_binary]].dropna()
        df_pred_gt[col_gt_binary] = df_pred_gt[col_gt_binary].astype(np.int8)
        df_pred_agg = pd.concat([df_pred_agg,df_pred_gt])

    #
    # Calculate ROC curve for given antibody, based on different thresholds on the continuous predictions compared to binary ground truth values
    col_pred = "{} {}".format(col_label_predictions, antibody)

    # Return auc and threshold table to be plotted later
    print("Calculating ROC curve and AUC (area under the curve)")
    return threshold(df_pred_agg[col_gt_binary], df_pred_agg[col_pred], scores)


def evaluate_elisa(path_data, path_elisa, protein, antibody, dilution_antibody):

    col_label_predictions = "Ratio positive"

    # Load ELISA data
    print('Loading ELISA data from {}'.format(path_elisa))
    df_elisa = pd.read_csv(path_elisa)
    df_elisa.rename(columns = {'sample name':'Sample Name', 'sample':'Sample'}, inplace = True) 
    elisa_cols = [c for c in df_elisa.columns if (('IgM' in c or 'IgG' in c or 'IgA' in c) and protein in c)]
    df_elisa = df_elisa[['Sample Name','Sample'] + elisa_cols]

    print('Loading per-well predictions and comparing to ELISA')
    df_pred_agg = pd.DataFrame()    
    file_list = [fn for fn in os.listdir(Path(path_data)) if (fn.endswith("csv") and (antibody in fn))]
    for filename in file_list:
        file_path = Path(path_data) / filename

        # Load per-well predictions
        df_predictions = pd.read_csv(file_path)
              
        # Exclude wells that have not transferred
        wells_excluded = (df_predictions.Transferred==False).sum()
        if wells_excluded>0:
            print("Excluding {} wells that have not transferred from evaluation".format(wells_excluded))
        df_predictions = df_predictions[df_predictions.Transferred==True]
        pred_cols = ["well", "Destination Plate Barcode", "Sample Name", "Sample", "Dilution", "Transferred"] + [c for c in df_predictions.columns if ('positive' in c)]
        df_predictions = df_predictions[pred_cols]

        #
        # Calculate dataset for per-channel continuous predictions       
        col_pred = "{} {}".format(col_label_predictions, antibody)
        col_elisa = "{} {}".format(protein, antibody)
        col_elisa_raw = "Raw {} {}".format(protein, antibody)

        # Only retain per-well ratios for wells that contain a sample for which ground truth exists explicitly
        df_pred_ch = df_predictions[df_predictions.Dilution == dilution_antibody[antibody]]
        df_pred_elisa = df_pred_ch.merge(df_elisa, on=['Sample Name', "Sample"], how='inner')
        df_pred_elisa[col_elisa_raw].replace('', np.nan, inplace=True) # Makes sure to exclude any samples for which no elisa data exists
        df_pred_elisa = df_pred_elisa[['Sample Name', 'Sample', col_pred, col_elisa, col_elisa_raw]].dropna()
        df_pred_elisa[col_elisa] = df_pred_elisa[col_elisa].astype(np.int8)
        df_pred_agg = pd.concat([df_pred_agg,df_pred_elisa])

    #
    # Calculate correlations between raw ELISA data and positivity rate predictions
    col_pred = "{} {}".format(col_label_predictions, antibody)
    col_elisa_raw = "Raw {} {}".format(protein, antibody)
    data_cov = df_pred_agg[df_pred_agg['Sample Name'].str.contains('S2020')]
    data_uu = df_pred_agg[df_pred_agg['Sample Name'].str.contains('S2017')]

    # Average correlations over 1000 subsampling runs, picking same number of neg as pos samples
    c_spearman = np.zeros(1000)
    for i in range(0,1000):
        data_uu_reduced = data_uu.sample(n=data_cov.shape[0])
        data_corr = pd.concat([data_cov,data_uu_reduced], axis=0)
        c_spearman[i] = data_corr[[col_pred,col_elisa_raw]].corr(method="spearman").iloc[0,1]
    c_spearman_mean = c_spearman.mean()
    c_spearman_std = c_spearman.std()

    data = df_pred_agg[['Sample Name','Sample',col_pred,col_elisa_raw]].copy()
    data = data.sort_values(by=col_elisa_raw)

    return (c_spearman_mean, c_spearman_std, data)
