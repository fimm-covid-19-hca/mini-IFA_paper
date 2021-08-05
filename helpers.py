"""
Analysis code for article: Image-Based & Machine Learning-Guided Multiplexed Serology Test for SARS-CoV-2

Helper methods for quality control, normalisation and plotting

@author: Christian Guckelsberger, Lassi Paavolainen
"""
def select_normalisation_wells(normalisation_type, controls_pos, controls_neg):
    # Fit normalisers on control wells for dilution (i.e. labelled and unlabelled samples)
    if(normalisation_type == "PosNegBal"):
        print("Normalisation based on balancing positive/negative controls")
        controls_pos_transferred = controls_pos.loc[controls_pos.Transferred==True,'well']
        controls_neg_transferred = controls_neg.loc[controls_neg.Transferred==True,'well']
        pick_n_transferred_each = min(controls_pos_transferred.shape[0],controls_neg_transferred.shape[0])
        normalisation_controls_wells = controls_pos_transferred.sample(n=pick_n_transferred_each).tolist() + controls_neg_transferred.sample(n=pick_n_transferred_each).tolist()
    elif(normalisation_type == "NegOnly"):
        print("Normalisation based on negative controls only")
        normalisation_controls_wells = controls_neg.loc[controls_neg.Transferred==True,'well'].tolist()
    elif(normalisation_type == "PosNeg"):
        print("Normalisation based on all positive/negative controls")
        normalisation_controls_wells = controls_neg.loc[controls_neg.Transferred==True,'well'].tolist() + controls_pos.loc[controls_pos.Transferred==True,'well'].values.tolist()

    return normalisation_controls_wells


def well_qc(df_features, dilution):
    # Feature used for QC to check if controls are transferred correctly
    qc_feat = 'NUCLEUS.INTENSITY-MEAN.DAPI'
    # Std multiplier for QC feature. The range is mean +- qc_mult * std
    qc_mult = 1.0

    pd_mean = df_features.loc[df_features.Dilution==dilution,qc_feat].mean()
    pd_std = df_features.loc[df_features.Dilution==dilution,qc_feat].std()
    pdwell_mean = df_features.loc[df_features.Dilution==dilution,['well',qc_feat]].groupby('well').mean()
    flagged_wells = pdwell_mean[(pdwell_mean[qc_feat] < pd_mean - qc_mult*pd_std) | (pdwell_mean[qc_feat] > pd_mean + qc_mult*pd_std)].index.values

    return flagged_wells
    

def check_sufficient_transfer(controls_pos_transferred, controls_neg_transferred, flagged_wells):
    controls_pos_transferred.loc[controls_pos_transferred['well'].isin(flagged_wells),'Transferred'] = False
    controls_neg_transferred.loc[controls_neg_transferred['well'].isin(flagged_wells),'Transferred'] = False
    insufficient_controls = False
    if controls_pos_transferred.Transferred.sum() < controls_pos_transferred.shape[0]:
        if controls_pos_transferred.Transferred.sum() >=2:
            print("Only {} of {} positive control wells were transferred -- still sufficient. Continuing normalisation.".format(controls_pos_transferred.Transferred.sum(), controls_pos_transferred.shape[0]))
        else:
            print("Only {} of {} positive control wells were transferred -- too few.".format(controls_pos_transferred.Transferred.sum(), controls_pos_transferred.shape[0]))
            insufficient_controls = True

    if controls_neg_transferred.Transferred.sum() < controls_neg_transferred.shape[0]:
        if controls_neg_transferred.Transferred.sum() >=2:
            print("Only {} of {} negative control wells were transferred -- still sufficient. Continuing normalisation.".format(controls_neg_transferred.Transferred.sum(), controls_neg_transferred.shape[0]))
        else:
            print("Only {} of {} negative control wells were transferred -- too few.".format(controls_neg_transferred.Transferred.sum(), controls_neg_transferred.shape[0]))
            insufficient_controls = True

    return not insufficient_controls
