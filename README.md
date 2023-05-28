# FPredX

This repository contains the FPredX models for the prediction of excitation maximum, emission maximum, brightness and oligomeric state of fluorescent proteins.

Tam, C.; Zhang, K.Y.J. FPredX: Interpretable models for the prediction of spectral maxima, brightness, and oligomeric states of fluorescent proteins. Proteins 2021, doi:10.1002/prot.26270. https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.26270

![png1](https://user-images.githubusercontent.com/51283097/172741590-a2c15ed7-1868-4a7d-b8f1-5b0ec3dd28aa.PNG)
![png2](https://user-images.githubusercontent.com/51283097/172741615-03148704-3815-412e-b3b5-c7b2616f381b.PNG)
![png3](https://user-images.githubusercontent.com/51283097/172742088-c4d9017e-6fb6-4133-9bb0-2aac73f70e15.PNG)

# Prerequisites

Step1: Please install the following

`MAFFT 7.471` (install MAFFT: https://mafft.cbrc.jp/alignment/software/source.html)

# Installation

Step2: Clone FPredX by

`git clone https://github.com/johnnytam100/FPredX.git`

Step3: Install required packages 

```
pip install --no-cache-dir -r requirements.txt
```

Done!

# Prediction
Simply 

`cd ./FPredX`

`cp (path to your fluorescent proteins fasta) ./`

`python FPredX_predict.py -m /usr/local/bin/mafft -p example.fasta`


# Remarks

**(1) Insertions that cannot fit into the aligned FP sequences of the original FPredX dataset were removed to maintain the total number of features.**

Higher the number of insertions removed, lower the prediction accuracy.

Check if any insertions were removed from the above log records of MAFFT OR confirm with the multiple aligned fasta file `FPredX_mafft_predict.fasta` generated during prediction.


**(2) If position-residue pairs from the sequences submitted for prediction were not present in the aligned FP sequences of the original FPredX dataset, all zeros were assigned to the positions for prediction.**

Higher the number of position-residue pairs reported, lower the prediction accuracy.

Check for locations of the position-residue pairs from the multiple aligned fasta file `FPredX_mafft_predict.fasta` generated during prediction.
