"""
Analysis code for article: Image-Based & Machine Learning-Guided Multiplexed Serology Test for SARS-CoV-2

Methods for plotting

@author: Christian Guckelsberger, Lassi Paavolainen
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_confusion_matrix(cm, classes):
    import pretty_print_cm as pp
    df_cm = pd.DataFrame(cm, index = [c for c in classes.values()], columns = [c for c in classes.values()])
    df_cm.rename(columns={"small bright":"small\nbright", "double positive": "double\npositive", "His positive":"HIS positive", "double negative": "double\nnegative"}, inplace=True)
    df_cm.index=df_cm.columns
    pp.pretty_plot_confusion_matrix(df_cm, percentages=False, fmt='d')


def plot_correlation(data, protein, antibody):

    fig_corr = plt.figure(figsize=(6,6))
    ax_corr = fig_corr.gca()   #Get Current Axis
    ax_corr.cla()

    # Plot correlations
    data_cov = data[data['Sample Name'].str.contains('S2020')]
    data_uu = data[data['Sample Name'].str.contains('S2017')]
    data_uu_reduced = data_uu.sample(n=data_cov.shape[0])
        
    # Export raw data used for correlation plot
    col_pred = "Ratio positive {}".format(antibody)
    col_gt_raw = "Raw {} {}".format(protein, antibody)
    ax_corr.plot(data_cov[col_gt_raw], data_cov[col_pred], 'o', alpha=0.5, color="red", label='COVID-19')
    ax_corr.plot(data_uu_reduced[col_gt_raw], data_uu_reduced[col_pred], 'o', alpha=0.5, color="darkcyan", label='Negative')
    ax_corr.set_xlabel('Raw ELISA Values',fontsize=14)
    ax_corr.set_ylabel('Predicted Positive Ratio',fontsize=14)
    #ax_corr.set_title('Positives={}, Negatives={}, Total={}'.format(data_cov.shape[0], data_uu_reduced.shape[0], data_cov.shape[0]+data_uu_reduced.shape[0]),fontsize=14)
    ax_corr.legend()

    fig_corr.show()


def plot_ROC_curve(scores_val):

    fig_roc = plt.figure(figsize=(6,6))
    ax_roc = fig_roc.gca() #Get Current Axis
    ax_roc.cla()

    ax_roc.plot((1-scores_val['specificity']), scores_val['sensitivity'], color="green", alpha=0.7, linewidth=3)
    ax_roc.plot([0, 1], [0, 1], transform=ax_roc.transAxes, color='black', linestyle='--')

    ax_roc.spines['right'].set_visible(False)
    ax_roc.spines['top'].set_visible(False)
    ax_roc.set_xlim([0.0, 1.01])
    ax_roc.set_ylim([0.0, 1.01])
    ax_roc.xaxis.set_ticks(np.arange(0, 1.01, 0.25))
    ax_roc.yaxis.set_ticks(np.arange(0, 1.01, 0.25))
    ax_roc.set_xlabel('False Positive Rate (1-Specificity)', fontsize=14)
    ax_roc.set_ylabel('True Positive Rate (Sensitivity)', fontsize=14)
    #ax_roc.set_title('{0}; GT: {1} p, {2} n, {3} total; AUC: {4:.2f}'.format(ch, counts[1],counts[0],y_gt.shape[0], area))
    
    fig_roc.show()
