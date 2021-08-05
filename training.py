#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis code for article: Image-Based & Machine Learning-Guided Multiplexed Serology Test for SARS-CoV-2

Methods for performing training and cross-validation with different classifiers

@author: Christian Guckelsberger, Lassi Paavolainen
"""
import uuid
import os
import shutil
import numpy as np
import joblib
import importlib
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from pathlib import Path
from sklearn.metrics import multilabel_confusion_matrix, confusion_matrix, make_scorer, accuracy_score
from sklearn.model_selection import ParameterGrid, cross_validate
from sklearn.dummy import DummyClassifier
from sklearn.utils.class_weight import compute_sample_weight

# Module-level variables to be shared between mp_training processes
X_train_shared = None
y_train_shared = None
X_validate_shared = None
y_validate_shared = None


'''
Multiprocessing training and evaluation of arbitrary estimator

Parameters wrapped into args:
    (ix, g): index and parameter dictionary 
    temp_dir: directory to store temporary model instance
    X_train, X_validate, y_train, y_validate: training and validation data
    sample_weights_train: sample weights for training data
    scores: dictionary of score functions to evaluate the trained model on

Returns:
    values on individual scores, evaluated on validation set 
'''
def mp_training(args):
    (ix, g), estimator, param_fixed, temp_dir, sample_weights_train, scores = args[0], args[1][0], args[1][1], args[1][2], args[1][3], args[1][4]
    clf = estimator(**g)
    clf.set_params(**param_fixed)
    if estimator.__name__ != "MLPClassifier":
        clf.fit(X_train_shared, y_train_shared, sample_weight = sample_weights_train)
    else:
        clf.fit(X_train_shared, y_train_shared)

    y_pred = clf.predict(X_validate_shared)

    cm = confusion_matrix(y_validate_shared, y_pred)

    scores_val = dict()
    for score in scores:
        scores_val[score] = scores[score]._score_func(y_validate_shared, y_pred)
        if score != "accuracy" and np.unique(y_pred).shape[0]>2:
            score_individual = scores[score]._score_func(y_validate_shared, y_pred, average='individual')
            scores_val["{} pos".format(score)] = score_individual[0]

    filename_temp_model = temp_dir / '{0}.joblib'.format(ix)
    joblib.dump(clf, filename=filename_temp_model, compress=True)
    del(clf)
    return (scores_val, cm)


'''
Sensitivity classification score: The proportion of actual positives that are correctly identified as such.
They can also take true negatives, false positives, false negatives and true positives as input (default: None)
This allows to compute the confusion matrix outside only once for the calculation of many different metrics
'''
def score_sensitivity(y_true, y_pred, average=None, label_pos=1):

    if average is None:
        average = 'binary' if ((np.unique(y_validate_shared).shape[0] == 2) or (np.unique(y_pred).shape[0]==2)) else 'micro'

    #print("average: {}".format(average))
    if average == 'micro':
        MCM = multilabel_confusion_matrix(y_true, y_pred)
        tp = MCM[:, 1, 1].sum()
        fn = MCM[:, 1, 0].sum()
    elif average == 'individual':
        MCM = multilabel_confusion_matrix(y_true, y_pred)
        tp = MCM[:, 1, 1]
        fn = MCM[:, 1, 0]
    elif average == 'binary':
        MCM = multilabel_confusion_matrix(y_true, y_pred, labels=[label_pos])
        tp = MCM[0, 1, 1]
        fn = MCM[0, 1, 0]
        
    return tp / (tp+fn)


'''
Accuracy classification score
'''
def score_accuracy(y_true, y_pred, MCM=None, average=None):
    return accuracy_score(y_true, y_pred)


'''
Specificity classification score: The proportion of actual negatives that are correctly identified as such.
'''
def score_specificity(y_true, y_pred, average=None, label_pos=1):
    
    if average is None:
        average = 'binary' if ((np.unique(y_validate_shared).shape[0] == 2) or (np.unique(y_pred).shape[0]==2)) else 'micro'
        
    if average == 'micro':
        MCM = multilabel_confusion_matrix(y_true, y_pred)
        tn = MCM[:, 0, 0].sum()
        fp = MCM[:, 0, 1].sum()
    elif average == 'individual':
        MCM = multilabel_confusion_matrix(y_true, y_pred)
        tn = MCM[:, 0, 0]
        fp = MCM[:, 0, 1]
    elif average == 'binary':
        MCM = multilabel_confusion_matrix(y_true, y_pred, labels=[label_pos])
        tn = MCM[0, 0, 0]
        fp = MCM[0, 0, 1]

    return tn / (tn+fp)


'''
Specsens classification score: Composite measure of specificity and sensitivity which retains the original range between [0,1].
SpecSens = Specificity x Sensitivity
'''
def score_specsens(y_true, y_pred, average=None, label_pos=1):
    return (score_specificity(y_true, y_pred, average, label_pos) * score_sensitivity(y_true, y_pred, average, label_pos))


'''
Provides the scores to base the assessment/tuning of a certain classifier on
'''
def setup_scores():
    # Define scores to report on
    scores = {'specificity':make_scorer(score_specificity),
            'sensitivity':make_scorer(score_sensitivity),
            'specsens':make_scorer(score_specsens),
            'accuracy':make_scorer(score_accuracy)
            }

    # Define score for which the model should be tuned
    score_tuning = 'specsens'

    return (scores, score_tuning)
        

'''
    Provide fixed and grid hyperparameters for the respective classifier.
    If hyperparameters are already given, use those. Otherwise use hardcoded parameter grid.
'''
def setup_hyperparameters(classifier_type, hyperparameters=None):

    hp_grid = dict()
    if hyperparameters is not None:
        hp_grid = hyperparameters
    elif classifier_type == "RF":
        # Grid from original article:
        # n_estimators: The number of trees in the forest. Usually, the larger the better. Try to find cut-off for optimal training time.
        # max_features: The size of the random subsets of features to consider when splitting a node. The lower the greater the reduction of variance, but also the greater the increase in bias.
        hp_grid  = [{'n_estimators': [50, 100, 200, 300, 500], 
                     'max_features': [None, 'sqrt']}] # I.e. either n_features or sqrt(n_features)

        # Downsized grid:
        hp_grid  = [{'n_estimators': [50, 300], 
                     'max_features': [None, 'sqrt']}] # I.e. either n_features or sqrt(n_features)

    elif classifier_type == "SVM":
        
        # Grid from original article:
        # C: regularisation parameter. Positive, default = 1.
        # gamma: kernel coefficient only for rbf and poly
        hp_grid  = [{'C': [0.1, 1, 10, 50, 100], 
                     'gamma': [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]}]

        # Downsized grid:
        hp_grid  = [{'C': [1, 10, 50], 
                     'gamma': [0.0005, 0.001]}]
                
    elif classifier_type == "ANN":

        # Grid from original article:
        # hidden_layer_sizes: neural network architectures with one to three hidden layers of varying size
        # alpha: regularization parameter
        # hp_grid  = [{'hidden_layer_sizes':[(64,32), (32,32), (128), (64), (32), (256,128,64),(128,64,32),(256,128),(128,64)],
        #              'alpha': [1.0, 0.1, 0.05, 0.01, 0.005, 0.001]}]

        # Downsized grid:
        hp_grid  = [{'hidden_layer_sizes':[(32, 32), (64), (32), (128,64)],
                      'alpha': [1.0, 0.1, 0.05, 0.001]}]
        
    hp_fixed = dict()
    if classifier_type == "SVM":
        hp_fixed = {'kernel': 'rbf',
                    'cache_size': 3500, 
                    'tol':1e-3, 
                    'max_iter': 500000}
            
    elif classifier_type == "ANN":
        hp_fixed = {'activation':'logistic',
                    'max_iter': 1000000, 
                    'learning_rate_init': 0.0001}
        
    return (hp_fixed, hp_grid)


def get_classifier(classifier_type):

    if classifier_type == "RF":
        return getattr(importlib.import_module("sklearn.ensemble"), "RandomForestClassifier")
    elif classifier_type == "SVM":
        return getattr(importlib.import_module("sklearn.svm"), "SVC")
    elif classifier_type == "ANN":
        return getattr(importlib.import_module("sklearn.neural_network"), "MLPClassifier")


'''
Load cross-validation data stored in train/validate format
'''
def load_data_cv(path_data, antibody):
    
    print('Loading cross-validation data features and labels from {}'.format(path_data))
    data = dict()
    file_list = [fn for fn in os.listdir(Path(path_data)) if (fn.endswith("hdf5") and (antibody in fn))]
    for filename in file_list:
        filepath = Path(path_data) / filename
        
        # Read training/validation data from .hdf5 into Pandas DataFrame
        X_train = pd.read_hdf(filepath, key='X_train')
        X_validate = pd.read_hdf(filepath, key='X_validate')
        y_train = pd.read_hdf(filepath, key='y_train')
        y_validate = pd.read_hdf(filepath, key='y_validate')
        
        data[filename] = (X_train.to_numpy(), X_validate.to_numpy(), y_train.Label.to_numpy(), y_validate.Label.to_numpy(), X_train.columns)

    return data


'''
Load training data
'''
def load_data(path_data, antibody):
    
    print('Loading data features and labels')
    
    filename = [fn for fn in os.listdir(Path(path_data)) if (fn.endswith("hdf5") and (antibody in fn))][0]
    filepath = Path(path_data) / filename

    # Read training data from .hdf5 into Pandas data.frame
    X = pd.read_hdf(filepath, key='X').to_numpy()
    y = pd.read_hdf(filepath, key='y').Label.to_numpy()
        
    return X, y


'''
Single-threaded training of a single classifier
'''    
def train_regular(X, y, classifier_type, hp_fixed, hp_grid):

    # Calculate sample weight vector to compensate unbalanced labels in dataset
    sample_weights_train = compute_sample_weight('balanced', y)
    unique, counts = np.unique(y, return_counts=True)

    hp = dict(list(hp_fixed.items()) + list(hp_grid.items()))
    classifier = get_classifier(classifier_type)()
    classifier.set_params(**hp)

    print("Training classifier {} with hyperparameters {}".format(classifier_type, hp))

    # Train model          
    if classifier_type != "ANN":
        classifier.fit(X, y, sample_weight = sample_weights_train)
    else:
        classifier.fit(X, y)

    # Return model
    print("Training completed.")
    return classifier


'''
Train and tune binary classifier on parameter grid.
Trains and evaluates classifiers on separate sets in parallel, using as many cores as are available on the system.
Trained models are temporarily written to the hard-drive.
Parameters are tuned based on the hyperparameter grid. Scores are calculated for all items in scores. The best classifier is chosen based on score_tuning. 

Returns the best trained classifier based on evaluation on the validation set 
'''
def train_and_tune(classifier_type, X_train, X_validate, y_train, y_validate, hp_grid, hp_fixed, classes, scores, score_tuning):
    
    # Calculate sample weight vector to compensate unbalanced labels in dataset
    sample_weights_train = compute_sample_weight('balanced',y_train)
    unique, counts = np.unique(y_train, return_counts=True)
                            
    # Establish baseline by training and evaluating two types of dummy classifiers
    scores_dummy = dict()
    for st in ('uniform', 'stratified'):
        dummy_clf = DummyClassifier(strategy=st)
        dummy_clf.fit(X_train, y_train)
        y_pred = dummy_clf.predict(X_validate)
        scores_dummy[st] = dict()
        for score in scores:
            scores_dummy[st][score] = scores[score]._score_func(y_validate, y_pred)
            
    print("Tuning hyper-parameters for {}".format(score_tuning))
    pg = ParameterGrid(hp_grid)
    scores_param = dict()
    for score in scores:
        scores_param[score] = np.zeros(len(pg))
        if score != "accuracy" and np.unique(y_validate).shape[0]>2:
            scores_param["{} pos".format(score)] = np.zeros(len(pg))
    
    # Create temporary directory to store trained models so they won't have to be kept in memory
    temp_dir = Path(__file__).parent / 'temp_models_{}'.format(uuid.uuid1())
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    # Set training and valudation data as module-level variables, so they can be shared between parallel processes
    # No need to duplicate these potentially big arrays many times
    # Python shares all module-level variables between multiprocessing threads
    global X_train_shared, X_validate_shared, y_train_shared, y_validate_shared
    X_train_shared = X_train
    X_validate_shared = X_validate
    y_train_shared = y_train
    y_validate_shared = y_validate

    # 'None': run as most as mp.cpu_count() asynchronous processes
    pool = mp.Pool(processes=None)
    classifier = get_classifier(classifier_type)
    args = [classifier, hp_fixed, temp_dir, sample_weights_train, scores]
    new_iterable = ([(ix[0],g), args] for ix, g in np.ndenumerate(pg))

    results = pool.map(mp_training, new_iterable)    

    pool.close()
    pool.join()

    cms = np.zeros((len(results), len(classes), len(classes)))
    for i in range(0, len(results)):
        scores = results[i][0]
        cms[i] = results[i][1]
        for score in scores:
            scores_param[score][i] = scores[score]
    
    param_ix_best = dict()
    for score in scores_param.keys():
        param_ix_best[score] = np.argmax(scores_param[score])

    best_params = ParameterGrid(hp_grid)[param_ix_best[score_tuning]]    
    print("Best parameter set for {}:".format(score_tuning))
    print(best_params)
    
    print("\nScores on validation set:")
    for score in scores:
        print("{}:".format(score))
        
        param_scores = []
        for i, g in np.ndenumerate(pg):   
            param_scores.append((scores_param[score][i], g))
        sorted_params = sorted(param_scores, key=lambda tup: tup[0], reverse=True)

        for (s, p) in sorted_params:
            print("{:.3f} for {}".format(s, p))

    return scores_param, cms, temp_dir, param_ix_best


'''
Perform cross-validation

The parameter data is a list of Pandas Data.Frames, each corresponding to one cross-validation fold and each implementing a split into training and validation set.
Parameters are tuned based on the hyperparameter grid. Scores are calculated for all items in scores. The best classifier is chosen based on score_tuning. 

Returns the highest-scoring parameter combination, and the corresponding confusion matrix, averaged over all folds.
'''
def cross_validate(data, classifier_type, hp_fixed, hp_grid, classes, scores, score_tuning):
            
    temp_dirs = dict()
    df_scores_fold = pd.DataFrame(columns=[*scores] + ["fold","hyperparameter index"])
    mcms_agg = None
    for file, d in data.items():
        X_train, X_validate, y_train, y_validate, _ = d
        print("\nTraining and validating {} on fold {}".format(classifier_type, file))
        print("Train/validation split: {} / {}".format(X_train.shape[0], X_validate.shape[0]))

        # Train and tune;
        # Store scores over all classifiers, and confusion matrices 
        # Not delete any models yet 
        scores_param, mcms, temp_dir, _ = train_and_tune(classifier_type, X_train, X_validate, y_train, y_validate, hp_grid, hp_fixed, classes, scores, score_tuning)

        # Keep record of all temporary directories with model files
        temp_dirs[file] = temp_dir

        # Record scores for this fold
        df_scores = pd.DataFrame(scores_param)
        df_scores["fold"] = "{}".format(os.path.basename(file))
        df_scores["hyperparameter index"] = np.arange(len(list(ParameterGrid(hp_grid))))
        df_scores_fold = df_scores_fold.append(df_scores,ignore_index=True)
        if mcms_agg is None:
            mcms_agg = mcms
        else:
            mcms_agg+=mcms

    # Determine best architecture for given score 
    df_scores_fold_mean = df_scores_fold.groupby("hyperparameter index",as_index=False).mean()
    param_ix_best = df_scores_fold_mean.at[np.argmax(df_scores_fold_mean[score_tuning]),"hyperparameter index"]
    param_best = [*ParameterGrid(hp_grid)][param_ix_best]
    print("\nModel with best mean scores on {}: {}".format(score_tuning, param_best))

    hp = []
    for x in df_scores_fold["hyperparameter index"]:
        hp.append(list(ParameterGrid(hp_grid))[x])
    df_scores_fold["hyperparameters"] = hp

    # Fetch summed confusion matrix for the best model
    cm = mcms_agg[param_ix_best].astype(np.int32)

    # Delete models / clear temp dirs
    for d in temp_dirs.values():
        shutil.rmtree(d)

    return param_best, cm
  

'''
Train a classifier on a single data set consisting of features and labels.
Requires the classifier type as string, a path to the data, a classes and a hyperparameter dictionary.
Returns the trained classifier.
'''
def train(classifier_type, path_data, antibody, hyperparameters):
    hp_fixed, hp_grid = setup_hyperparameters(classifier_type, hyperparameters)
    X, y = load_data(path_data, antibody) 
    return train_regular(X, y, classifier_type, hp_fixed, hp_grid)


'''
Train and validate a classifier on a data set with a train/validation for different hyperparameters.
Requires the classifier type as string, a path to the data, and a classes dictionary.
Returns the highest-scoring parameter combination, and the corresponding confusion matrix, averaged over all folds.
'''
def train_cv(classifier_type, path_data, antibody, classes):
    scores, score_tuning = setup_scores()
    hp_fixed, hp_grid = setup_hyperparameters(classifier_type)
    data = load_data_cv(path_data, antibody) 
    return cross_validate(data, classifier_type, hp_fixed, hp_grid, classes, scores, score_tuning)
