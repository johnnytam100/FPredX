# FPredX

This repository contains the FPredX models for the prediction of excitation maximum, emission maximum, brightness and oligomeric state of fluorescent proteins.

# Prerequisites

Step1: Please install the following

`MAFFT 7.471` (install MAFFT: https://mafft.cbrc.jp/alignment/software/source.html)

`pip install joblib==1.1.0`

`pip install dill==0.3.4`

`pip install xgboost==0.90`

`pip install scikit-learn==0.22.2`

`pip install pandas==1.1.5`

`pip install numpy==1.21.2`

# Installation

Step2: clone FPredX by

`git clone https://github.com/johnnytam100/FPredX.git`

Step3: Change the path inside `FPredX_predict.py` to your local MAFFT executable.

Done!

# Prediction
Simply 

`cd ./FPredX`

`cp (path to your fluorescent proteins fasta) ./`

`python FPredX_predict.py (your fluorescent proteins fasta)`

# Usage example

`cd ./FPredX`

`python FPredX_predict.py example.fasta`

# Remarks

**(1) Insertions that cannot fit into the aligned FP sequences of the original FPredX dataset were removed to maintain the total number of features.**

Higher the number of insertions removed, lower the prediction accuracy.

Check if any insertions were removed from the above log records of MAFFT OR confirm with the multiple aligned fasta file `FPredX_mafft_predict.fasta`


**(2) If position-residue pairs from the sequences submitted for prediction were not present in the aligned FP sequences of the original FPredX dataset, all zeros were assigned to the positions for prediction.**

Higher the number of position-residue pairs reported, lower the prediction accuracy.

Check for locations of the position-residue pairs from the multiple aligned fasta file `FPredX_mafft_predict.fasta`
