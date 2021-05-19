# Scripts used structuraly analyze breast carcinogen chemicals
Project developed to analysis chemicals that can be classified as a breast carcinogens or related to 

# Dependancy
## Coding language 
- Python3.9
- Python library => CompDesc, bitarray, RDKIT
- R 4.0
- RDKIT (2020) in environement conda rdkit-env

## Prediction - QSPR models 
- genotoxicity: TEST (https://www.epa.gov/chemical-research/toxicity-estimation-software-tool-test#install)
- OPERA ER-AR prediction

# To do list
## Priority 1
- make clustering using toxifor on BCC and P4up/E2up chemicals
- cluster using chemicals type classification 
- cluster using the toxiprint
- reproduce clustering with the 10 hormone tested in kaumas
- Add pca for the chemicals overlap in addition to the upset plot
- ~~do a clustering with the 10 hormones from Agnes paper~~
- check link with SP data in the dropbox folder
- ~~Add borderline in the E2 up and P4 up !!!~~
- change CR value by the file - H295R_AllChemsTested_ANOVA_directionality_filtered_output_2018-08-09.xlsx
- ~~change hormone order OHPREG - PROG - OHPREG - DOC - CORTISOL - 11DC .. from kaumas 2016~~
- AFC with E2 and p4 - take a look on complete matrix

## Priority 2
- link to the dropbox SP with the data input
- check metric for the toxprint
- use a fingerprint based on the structural FP

# Main updates to report
- 05-03-21: init project from chemicals clustering - clean up legacy functions 
