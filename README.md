# mini-IFA_paper

We provide, notebook.ipynb a self-contained Jupyter notebook script which comprises the full ML model training, prediction and evaluation pipeline presented in the manuscript. It allows for the reproduction of our findings, based on the configuration in the "Settings" cell. Reproducing all reported findings requires multiple runs of the notebook with different settings, but it is configured by default to reproduce our most significant findings as reported in the manuscript.

The notebook utilizes various custom Python modules which are included in the root directory. It moreover makes use of metadata files and raw data, which are stored in the data_auxiliary and data_raw directories, respectively.

## System requirements

The Python code has been implemented and tested on Ubuntu 18.04 Linux with Python 3.6.9. Our virtual development server features a 16 core AMD EPYC processor (2.0GHz) with 64GB RAM. The code should work with any operating system with Python 3.6.x or newer that has sufficient RAM to handle the image feature data in parallel training processes. We confirm that the code works with the following versions of the following Python packages:
- joblib 0.14.1
- jupyter 1.0.0
- matplotlib 3.2.1
- numpy 1.18.4
- pandas 1.1.5
- scikit-learn 0.24.1
- seaborn 0.10.1

The notebook has been tested on the following Jupyter components:
- jupyter core 4.7.1
- jupyter-notebook 6.4.0
- jupyter client   : 6.1.12

Based on a conservative estimate, one pass through the notebook takes about 60 minutes, given the hardware specifications and library versions as indicated above.

## Quality Control (QC) and Visualisation

The QCandVis directory comprises the R scripts used to run the QC pipeline and to create all result figures as shown in the paper. Most figures can also be obtained directly from the Jupyter notebook, but with a slightly different optic. Follow the steps below to perform the QC and data visualisation:

1. Create the following directory structure and set "/Results/Normalisation_NegOnly" as the working directory in R:
- /Data/Images
- /Data/Normalisation_NegOnly           
- /Results/Normalisation_NegOnly
- /Code
- /RDA

2. Run the following scripts in sequence:
- Misc_Functions.R
- Data_processing.R
- QC_Plots.R
- Images_Scatterplots.R
- Analysis_Plots.R

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
