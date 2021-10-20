from argparse import ArgumentParser
import os
import numpy as np
import pandas as pd
import glob
import joblib

parser = ArgumentParser()
parser.add_argument('arg1')
args = parser.parse_args()

###################################
''' Multiple Sequence Alignment '''
###################################

mafft_path = "/bwfefs/home/cltam/script/mafft-build/bin/mafft"

os.system(mafft_path + " --add " + args.arg1 + " --keeplength FPredX_mafft.fasta > FPredX_mafft_predict.fasta")

##########################
''' Fasta to dataframe '''
##########################

### Load fasta
filein = open('FPredX_mafft_predict.fasta', 'r', encoding='utf-8')
fasta = filein.readlines()
filein.close()

### Convert multiple lines fasta to single line fasta
for i in range(len(fasta)-1):
  if fasta[i].count('>') == 0 and fasta[i+1].count('>') == 0:
    fasta[i] = fasta[i].replace('\n','')

for i in range(len(fasta)-2, -1, -1):
  if fasta[i][-1:] != "\n":
    fasta[i] = fasta[i] + fasta.pop(i+1)

### Create an empty dataframe to collect seq name and sequences
seq_list_df = pd.DataFrame()

### Prepare seq_name list 
name = []
for i in range(len(fasta[::2])):
  name.append(fasta[::2][i].replace('\n', '').replace('>',''))
seq_list_df['name'] = name

### Prepare sequence list
seq_align = []
for i in range(len(fasta[1::2])):
  seq_align.append(fasta[1::2][i].replace('\n', ''))
seq_list_df['seq_align'] = seq_align


########################
''' One hot encoding '''
########################

byposition_pred = seq_list_df['seq_align'].apply(lambda x:pd.Series(list(x)))
one_hot_df_pred = pd.get_dummies(byposition_pred)

### Remove unavailable residue
available_res = pd.read_csv("available_res.csv", index_col=0)
one_hot_df_pred_trim = pd.DataFrame()
unavailable_list = []
for i in one_hot_df_pred.columns:
  if i in available_res.index:
    one_hot_df_pred_trim = pd.concat([one_hot_df_pred_trim, one_hot_df_pred[i]], axis=1)
  else:
    unavailable_list.append(i)

### Remove FPredX seq
num_pred = int(os.popen("grep -c '>' " + args.arg1).read())
seq_list_df_trim = seq_list_df.iloc[-num_pred:,:]
one_hot_df_pred_trim = one_hot_df_pred_trim.iloc[-num_pred:,:]

one_hot_df_pred_trim.reset_index(drop = True, inplace = True)
seq_list_df_trim.reset_index(drop = True, inplace = True)

###############################
''' Predict with all models '''
###############################

all_pred_mean_df = pd.DataFrame()

dir_list = ["ex_model",
            "em_model",
            "bright_model",
            "oligo_model"]

# Results for append to final table
mut_pred_df_list = []

# Results for mut-wt difference boxplot plotting
global mut_pred_df_all


for dir in dir_list:

  if dir != "oligo_model":

    mut_pred_df_all = pd.DataFrame()

    model_list = glob.glob("./" + dir + "/*")

    model_list.sort()

    for model in model_list:

      # predict
      xgb = joblib.load(model)

      # Results for append to final table
      mut_features = np.array(one_hot_df_pred_trim)
      mut_pred = xgb.predict(mut_features)
      mut_pred_df = pd.DataFrame(mut_pred)
      mut_pred_df.columns = [(dir + model[-2:] + "_pred")]
      mut_pred_df_all = pd.concat([mut_pred_df_all, mut_pred_df], axis=1)

    # calculate mean
    mean = pd.DataFrame(mut_pred_df_all.mean(axis=1), columns=[dir])

    # concat, save
    all_pred_mean_df = pd.concat([all_pred_mean_df, mean], axis=1)

    mut_pred_df_all_merge = pd.concat([mut_pred_df_all, seq_list_df_trim, one_hot_df_pred_trim], axis=1)
    mut_pred_df_all_merge.to_csv(dir + "_pred.csv")

  else:

    mut_pred_df_all = pd.DataFrame()

    model_list = glob.glob("./" + dir + "/*")

    model_list.sort()
    
    for model in model_list:

      # predict
      xgb = joblib.load(model)

      # Results for append to final table
      mut_features = np.array(one_hot_df_pred_trim)
      mut_pred = xgb.predict(mut_features)
      mut_pred_df = pd.DataFrame(mut_pred)
      mut_pred_df.columns = [(dir + model[-2:] + "_pred")]
      mut_pred_df_all = pd.concat([mut_pred_df_all, mut_pred_df], axis=1)

    # calculate mean
    mean = mut_pred_df_all.mode(axis=1)[0]
    mean.rename(dir, inplace=True)

    # concat, save
    all_pred_mean_df = pd.concat([all_pred_mean_df, mean], axis=1)

    mut_pred_df_all_merge = pd.concat([mut_pred_df_all, seq_list_df_trim, one_hot_df_pred_trim], axis=1)
    mut_pred_df_all_merge.to_csv(dir + "_pred.csv")

### Save prediction summary
all_pred_mean_df = pd.concat([all_pred_mean_df,  seq_list_df_trim, one_hot_df_pred_trim], axis=1)
all_pred_mean_df.to_csv("summary_pred_mean.csv")

### Report unavailable position-residue pairs

if len(unavailable_list) != 0:
  print("\n\n#####################")
  print(">>> FPredX NOTICE <<<")
  print("#####################")
  print("\n\n(1) Insertions that cannot fit into the aligned FP sequences of the original FPredX dataset were removed to maintain the total number of features.")
  print("\nHigher the number of insertions removed, lower the prediction accuracy.")
  print("\nCheck if any insertions were removed from the above log records of MAFFT OR confirm with the multiple aligned fasta file 'FPredX_mafft_predict.fasta'")
  print("\n\n(2) Position-residue pairs:\n")
  print(unavailable_list)
  print("\nfrom the sequences submitted for prediction were not present in the aligned FP sequences of the original FPredX dataset.")
  print("\nAll zeros were assigned to the positions for prediction.")
  print("\nHigher the number of position-residue pairs reported above, lower the prediction accuracy.")
  print("\nCheck for locations of the position-residue pairs from the multiple aligned fasta file 'FPredX_mafft_predict.fasta'")
  print("\n\nPrediction complete. Thanks for using FPredX! ( ^_____^  )\n\n")

else:
  print("\n\n#####################")
  print(">>> FPredX NOTICE <<<")
  print("#####################")
  print("\n\n(1) Insertions that cannot fit into the aligned FP sequences of the original FPredX dataset were removed to maintain the total number of features.")
  print("\nHigher the number of insertions removed, lower the prediction accuracy.")
  print("\nCheck if any insertions were removed from the above log records of MAFFT OR confirm with the multiple aligned fasta file 'FPredX_mafft_predict.fasta'")
  print("\n\nPrediction complete. Thanks for using FPredX! ( ^_____^  )\n\n")
