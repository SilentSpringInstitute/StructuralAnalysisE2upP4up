# Scripts used to build E2up and P4up QSAR models and to use them to predict E2up and P4up properties on a set mammary carcinogens 
All scripts are included in the folder sources and are separated by type of scripts 
- R scripts, mostly used to draw figures
- Python scripts used for data analysis


# Dependencies
## By coding languages 
- Python3.9
    - RDKIT (2020) in environement conda rdkit-env
    - CompDesc ($pip install -i https://test.pypi.org/simple/ CompDesc)
    - scipy ($pip install scipy)
    - molvs ($pip install molvs)
    - bitarray ($pip install bitarray)
    - openpyxl ($pip install openpyxl)
    - tensorflow ($pip install tensorflow) use keras from tensorflow package
    - sklearn ($pip install sklearn)

- R 4.0
    - Toolbox from github (https://github.com/ABorrel/Toolbox_packageR)
    - Clone QSAR R scripts available at https://github.com/ABorrel/QSAR-QSPR

- set up the path for a internal run of the R script in the python script source/py/runExternal.py

## Prediction - QSAR models 
- genotoxicity: TEST (https://www.epa.gov/chemical-research/toxicity-estimation-software-tool-test#install)
- ER-AR prediction: OPERA2.7 (https://github.com/kmansouri/OPERA)
<br><br><br>

# Input data
The input folder included input file to reproduce the analysis:
- list of exposure "BCRelExposureSources_P65_051221.csv"
- list of chemicals considered "cross_lists_for_analysis_042522.xlsx"
- list of hormones "hormones.csv"
- QSAR models "QSARs/"
<br><br><br>

# Run the script
A Jupiter notebook was developed to run all analysis partially or completely and apply QSAR models. <br> 
[Jupiter notebook](https://github.com/SilentSpringInstitute/StructuralAnalysisE2upP4up/blob/master/sources/py/main.ipynb) 

<br><br><br>

