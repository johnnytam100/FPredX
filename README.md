# FPredX

This repository contains the FPredX models for the prediction of excitation maximum, emission maximum, brightness and oligomeric state of fluorescent proteins.

# Prerequisites

# Installation
No installation is required.
But change the path to your local MAFFT executable inside `FPredX_predict.py`.

# Prediction
Simply 

`python FPredX_predict.py (your fluorescent proteins fasta)`

# Usage example

`cd example`

`python FPredX_predict.py pred.fasta`

# Remarks

**(1) Insertions that cannot fit into the aligned FP sequences of the original FPredX dataset were removed to maintain the total number of features.**

Higher the number of insertions removed, lower the prediction accuracy.

Check if any insertions were removed from the above log records of MAFFT OR confirm with the multiple aligned fasta file `FPredX_mafft_predict.fasta`

**(2) If position-residue pairs from the sequences submitted for prediction were not present in the aligned FP sequences of the original FPredX dataset, all zeros were assigned to the positions for prediction.**

Higher the number of position-residue pairs reported above, lower the prediction accuracy.

Check for locations of the position-residue pairs from the multiple aligned fasta file `FPredX_mafft_predict.fasta`
