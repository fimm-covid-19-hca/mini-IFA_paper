#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis code for article: Image-Based & Machine Learning-Guided Multiplexed Serology Test for SARS-CoV-2

Fixed settings to be used for the correct operation of Jupiter Python Notebook.

@author: Christian Guckelsberger, Lassi Paavolainen
"""

from pathlib import Path

# Input and output data path
path_data_raw = Path("data_raw/")
path_data_auxiliary = Path("data_auxiliary/")
path_output = Path("data_out/")

# Filenames for classes, controls, metadata and original hyperparameters csv files (in input data folder)
file_classes = "classes.csv"
file_controls = "controls.csv"
file_metadata = "metadata.csv"
file_original_hyperparameters = "default_hyperparameters.csv"

# Filenames for evaluation comparison data (in input data folder)
file_elisa = "elisa.csv"
file_visual_gt = "visual_gt.csv"

# Dilutions to filter data from per antibody 
# Specified antibodies show stronger signals under specific dilutions, hence we focus on those dilutions for training and prediction 
dilution_antibody = {'IgM': 25, 'IgG': 100, 'IgA': 25}
