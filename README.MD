# Scripts used to build E2up and P4up QSAR models and to use them to predict E2up and P4up properties on a set mammary carcinogens 
All scripts are included in the folder sources and are separated by type of script 
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

## Prediction - QSPR models 
- genotoxicity: TEST (https://www.epa.gov/chemical-research/toxicity-estimation-software-tool-test#install)
- ER-AR prediction: OPERA2.7 (https://github.com/kmansouri/OPERA)
<br><br><br>

# Input data
Here put link to requiered input files
<br><br><br>

# Run the script
A Jupiter notebook was developed to run all analysis partially or completely. <br>
> sources/py/main.ipynb

<br><br><br>

# To do list
## Priority 1
- ~~make clustering using toxifor on BCC and P4up/E2up chemicals~~
- ~~cluster using chemicals type classification~~ 
- ~~cluster using the toxiprint~~
- ~~reproduce clustering with the 10 hormone tested in kaumas~~
- ~~Add pca for the chemicals overlap in addition to the upset plot~~
- ~~do a clustering with the 10 hormones from Agnes paper~~
- ~~check link with SP data in the dropbox folder~~
- ~~Add borderline in the E2 up and P4 up !!!~~
- ~~change CR value by the file H295R_AllChemsTested_ANOVA_directionality_filtered_output_2018-08-09.xlsx~~
- ~~change hormone order OHPREG - PROG - OHPREG - DOC - CORTISOL - 11DC .. from kaumas 2016~~
- ~~AFC with E2 and p4 - take a look on complete matrix~~
- ~~Need to change CR value in stereogensis with hagren2019~~
- ~~test Ghost approach for undersampling~~
- ~~add DNN for the classification prediction~~
- ~~remove borderline chemicals~~
- ~~project on the SOM different color~~ 
- ~~make a undersampling variable~~

## Priority 2
- ~~check metric for the toxprint~~
- ~~use a fingerprint based on the structural FP~~
- ~~Add combination of ToxPrint to improve building unit~~
