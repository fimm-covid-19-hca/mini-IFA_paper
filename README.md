# mini-IFA_paper

notebook.ipynb is a self-explanatory Jupyter notebook script presenting the full ML model training, prediction and evaluation pipeline used in the manuscript. It utilizes various custom Python modules included in the root directory. In addition, metadata files used in the work are included in the data_auxiliary directory.

## System requirements

Python code is implemented and tested in Ubuntu 18.04 with Python 3.6.9. Our virtual development server includes 16 core AMD EPYC processor (2.0GHz) with 64GB RAM. The code should work with any OS with Python 3.6.x or newer from 3.x series having enough RAM to handle the parallel training process. The code works at least with the following versions of Python packages:
- joblib 0.14.1
- jupyter 1.0.0
- matplotlib 3.2.1
- numpy 1.18.4
- pandas 1.1.5
- scikit-learn 0.24.1
- seaborn 0.10.1

## QC and Visualisation

QCandVis directory includes R scripts used to run QC pipeline and create the result figures.

Create following directory structure and set "/Results/Normalisation_NegOnly" as working directory:
- /Data/Images
- /Data/Normalisation_NegOnly           
- /Results/Normalisation_NegOnly
- /Code
- /RDA

Scripts to be run in the sequence to produce QC and data visualisation:
- Misc_Functions.R
- Data_processing.R
- QC_Plots.R
- Images_Scatterplots.R
- Analysis_Plots.R
