import pandas as pd
import numpy as np
import argparse
import glob
import joblib
import sys
import os

def parse(args=None):
    '''arguments.
    '''
    parser = argparse.ArgumentParser(prog='FPredX-Fluorescence prediction',
                                     description='Run FPredX models for the prediction of excitation maximum, emission maximum, brightness and oligomeric state of fluorescent protein filter-Agmata',
                                     epilog='(c) authors')
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s ' + '1',
                        help="Show program's version number and exit.")
    parser.add_argument('-m', '--mafft',
                        help='path to mafft installation',
                        required=True)
    parser.add_argument('-p', '--path',
                        help='path to fasta file',
                        required=True)
    parser.add_argument('-o', '--out',
                        help='path to output directory',
                        required=True)
    results = parser.parse_args(args)
    return results

def load_fasta(filename):
    filein = open(filename, 'r', encoding='utf-8')
    fasta = filein.readlines()
    filein.close()
    return fasta

def convert_multiline_fasta_to_singleline(fasta):
    for i in range(len(fasta)-1):
        if fasta[i].count('>') == 0 and fasta[i+1].count('>') == 0:
            fasta[i] = fasta[i].replace('\n','')

    for i in range(len(fasta)-2, -1, -1):
        if fasta[i][-1:] != "\n":
            fasta[i] = fasta[i] + fasta.pop(i+1)
    return fasta

def prepare_sequence_dataframe(fasta):
    seq_list_df = pd.DataFrame()
    name = []
    for i in range(len(fasta[::2])):
        name.append(fasta[::2][i].replace('\n', '').replace('>',''))
    seq_list_df['name'] = name
    seq_align = []
    for i in range(len(fasta[1::2])):
        seq_align.append(fasta[1::2][i].replace('\n', ''))
    seq_list_df['seq_align'] = seq_align
    return seq_list_df

def one_hot_encoding(seq_align_df):
    byposition_pred = seq_align_df['seq_align'].apply(lambda x:pd.Series(list(x)))
    one_hot_df_pred = pd.get_dummies(byposition_pred)
    return one_hot_df_pred

def remove_unavailable_residues(one_hot_df_pred):
    available_res = pd.read_csv("available_res.csv", index_col=0)
    one_hot_df_pred_trim = pd.DataFrame()
    unavailable_list = []
    for i in one_hot_df_pred.columns:
        if i in available_res.index:
            one_hot_df_pred_trim = pd.concat([one_hot_df_pred_trim, one_hot_df_pred[i]], axis=1)
        else:
            unavailable_list.append(i)
    return one_hot_df_pred_trim, unavailable_list

def remove_fpredx_sequences(one_hot_df_pred_trim, seq_list_df, num_pred):
    seq_list_df_trim = seq_list_df.iloc[-num_pred:,:]
    one_hot_df_pred_trim = one_hot_df_pred_trim.iloc[-num_pred:,:]
    one_hot_df_pred_trim.reset_index(drop = True, inplace = True)
    seq_list_df_trim.reset_index(drop = True, inplace = True)
    return seq_list_df_trim, one_hot_df_pred_trim

def predict_with_all_models(seq_list_df_trim, one_hot_df_pred_trim):
    all_pred_mean_df = pd.DataFrame()
    dir_list = ["ex_model", "em_model", "bright_model", "oligo_model"]
    mut_pred_df_list = []
    global mut_pred_df_all
    for dir in dir_list:
        if dir != "oligo_model":
            mut_pred_df_all = pd.DataFrame()
            model_list = glob.glob("./" + dir + "/*")
            model_list.sort()
            for model in model_list:
                xgb = joblib.load(model)
                mut_features = np.array(one_hot_df_pred_trim)
                mut_pred = xgb.predict(mut_features)
                mut_pred_df = pd.DataFrame(mut_pred)
                mut_pred_df.columns = [(dir + model[-2:] + "_pred")]
                mut_pred_df_all = pd.concat([mut_pred_df_all, mut_pred_df], axis=1)
            mean = pd.DataFrame(mut_pred_df_all.mean(axis=1), columns=[dir])
            all_pred_mean_df = pd.concat([all_pred_mean_df, mean], axis=1)
            mut_pred_df_all_merge = pd.concat([mut_pred_df_all, seq_list_df_trim, one_hot_df_pred_trim], axis=1)
            mut_pred_df_all_merge.to_csv(dir + "_pred.csv")
        else:
            mut_pred_df_all = pd.DataFrame()
            model_list = glob.glob("./" + dir + "/*")
            model_list.sort()
            for model in model_list:
                xgb = joblib.load(model)
                mut_features = np.array(one_hot_df_pred_trim)
                mut_pred = xgb.predict(mut_features)
                mut_pred_df = pd.DataFrame(mut_pred)
                mut_pred_df.columns = [(dir + model[-2:] + "_pred")]
                mut_pred_df_all = pd.concat([mut_pred_df_all, mut_pred_df], axis=1)
            mean = mut_pred_df_all.mode(axis=1)[0]
            mean.rename(dir, inplace=True)
            all_pred_mean_df = pd.concat([all_pred_mean_df, mean], axis=1)
            mut_pred_df_all_merge = pd.concat([mut_pred_df_all, seq_list_df_trim, one_hot_df_pred_trim], axis=1)
            mut_pred_df_all_merge.to_csv(dir + "_pred.csv")
    return all_pred_mean_df

def save_prediction_summary(all_pred_mean_df, seq_list_df_trim, one_hot_df_pred_trim, out_path):
    all_pred_mean_df = pd.concat([all_pred_mean_df, seq_list_df_trim, one_hot_df_pred_trim], axis=1)
    all_pred_mean_df.to_csv(out_path + "/summary_pred_mean.csv")
    
def report_unavailable_residue_pairs(unavailable_list):
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
        print("\n\nPrediction complete. Thanks for using FPredX! ( ^^ )\n\n")
    else:
        print("\n\n#####################")
        print(">>> FPredX NOTICE <<<")
        print("#####################")
        print("\n\n(1) Insertions that cannot fit into the aligned FP sequences of the original FPredX dataset were removed to maintain the total number of features.")
        print("\nHigher the number of insertions removed, lower the prediction accuracy.")
        print("\nCheck if any insertions were removed from the above log records of MAFFT OR confirm with the multiple aligned fasta file 'FPredX_mafft_predict.fasta'")
        print("\n\nPrediction complete. Thanks for using FPredX! ( ^^ )\n\n")
        
def run_mafft(mafft_path, fasta_path):
    os.system(mafft_path + " --add " + fasta_path + " --keeplength FPredX_mafft.fasta > FPredX_mafft_predict.fasta")

def main():
    
    args = parse(sys.argv[1:])
    print("Analysis Starts....")
    run_mafft(args.mafft, args.path)
    fasta = load_fasta('FPredX_mafft_predict.fasta')
    fasta = convert_multiline_fasta_to_singleline(fasta)
    seq_list_df = prepare_sequence_dataframe(fasta)
    one_hot_df_pred = one_hot_encoding(seq_list_df)
    one_hot_df_pred_trim, unavailable_list = remove_unavailable_residues(one_hot_df_pred)
    num_pred = int(os.popen("grep -c '>' " + args.path).read())
    seq_list_df_trim, one_hot_df_pred_trim = remove_fpredx_sequences(one_hot_df_pred_trim, seq_list_df, num_pred)
    all_pred_mean_df = predict_with_all_models(seq_list_df_trim, one_hot_df_pred_trim)
    save_prediction_summary(all_pred_mean_df, seq_list_df_trim, one_hot_df_pred_trim, args.out)
    report_unavailable_residue_pairs(unavailable_list)
    print('Analysis Finished!')

if __name__ == '__main__':
    main()
