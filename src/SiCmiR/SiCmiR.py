# coding=utf-8
import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import zscore
import torch
import torch.nn as nn

torch.set_default_dtype(torch.float)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#----Part I: Data preparation----#
def build_submatrix(df_all: pd.DataFrame, picked_names, picked_alias, output_dir, fill_value=0, save_extract=None):
    sub = pd.DataFrame(index=df_all.index)
    cols = []
    for i in range(977):
        name = picked_names[i]
        name_alias = picked_alias[i]
        if name in df_all.columns:
            col = df_all[name]
            col.name = name
        elif name_alias in df_all.columns:
            col = df_all[name_alias]
            col.name = name
        else:
            col = pd.Series(fill_value, index=df_all.index) #fill 0 if a required gene does not exist
            col.name = name   
        cols.append(col)
    sub = pd.concat(cols, axis=1)
    if save_extract:
        sub.to_csv(os.path.join(output_dir,save_extract))



    return sub

#----Part II: Prediction----#

class NeuralNet(nn.Module):
    def __init__(self, input_size=977, hidden_size=1024, output_size=1298):
        super(NeuralNet, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.ReLU(),
            nn.Dropout(p=0.3),  # dropout训
            nn.Linear(hidden_size, hidden_size),
                       
            nn.ReLU(),
            nn.BatchNorm1d(hidden_size),
            nn.Dropout(p=0.3),
            
            nn.Linear(hidden_size, output_size),
        )
    def forward(self, x):
        x = self.model(x)
        return x

def predict(
    model_path,
    input_path,
    output_path,
    output_dir='./',
    storage_dir='./data',
    extract=True,
    normalization=False,
    save_extract=None,
    save_zscore_input=None
):
    
    # Load model
    net = torch.load(f"{model_path}", map_location=device)
    net.eval()

    if extract:
        L1000_gene = pd.read_csv(os.path.join(storage_dir,'L1000gene.csv'), header=None)
        picked_names = L1000_gene.loc[:, 0]
        picked_alias = L1000_gene.loc[:, 1]
        data_matrix = pd.read_csv(input_path, index_col=0).T

        extracted_mRNA = build_submatrix(
            data_matrix, picked_names, picked_alias, output_dir, save_extract=save_extract)
    else:
        extracted_mRNA = pd.read_csv(input_path, index_col=0)

    if normalization:
        extracted_zscored_mRNA = extracted_mRNA.apply(lambda x: zscore(x, ddof=0), axis=0).fillna(0) #zscore scaling and fill NA with 0
        if save_zscore_input:
            extracted_zscored_mRNA.to_csv(os.path.join(output_dir,save_zscore_input))
        X_test_P = extracted_zscored_mRNA
    else:
        X_test_P = extracted_mRNA

    X_test_Tensor = torch.tensor(X_test_P.values, dtype=torch.float)
    miRNA_1298 = pd.read_csv(os.path.join(storage_dir,'1298miRNA.csv'), header=None).loc[:, 0]
    y_test_out = pd.DataFrame(net(X_test_Tensor).detach().numpy())
    y_test_out.columns = miRNA_1298
    y_test_out.index = X_test_P.index

    y_test_out.to_csv(os.path.join(output_dir,output_path))
    print(f"Results saved to {output_dir}")
    return y_test_out

def main():
    parser = argparse.ArgumentParser(description="SiCmiR: Predict miRNA expression from mRNA data")
    parser.add_argument("--input", "-i", type=str, default="./example/test_mRNA.csv",
        help="Input mRNA matrix (genes in rows, samples in columns). \
        Remove sample batch effects before use is preferable. Default: %(default)s")
    parser.add_argument("--output", "-o", type=str, default="predicted_miRNA.csv",
        help="Output predicted miRNA expression file. Default: %(default)s")
    parser.add_argument("--output_dir", "-od", default='./',
        help="Directory to save output files. Default: %(default)s")
    parser.add_argument("--storage_dir", "-sd", default='./data/',
        help="Directory to store pretrained model and auxiliary files. Default: %(default)s")
    parser.set_defaults(extract=True)
    parser.add_argument("--no_extract", "-n",dest="extract", action='store_false',
        help="Do not extract 977 landmark genes from input mRNA matrix (Not Recomended).")
    parser.add_argument("--normalization", "-norm", action='store_true',
        help="Whether to perform z-score normalization on input matrix. Default: False (No zscore normalization)")
    parser.add_argument("--save_extract", "-se", default=None,
        help="Output file name for extracted but unnormalized mRNA matrix (if `--extract` is used). Default: None (Not saving)")
    parser.add_argument("--save_zscore_input", "-z", default=None,
        help="Save zscored mRNA matrix used directly for prediction (if `--normalization` is used). Default: None (Not saving)")
    parser.add_argument("--model", "-m", default="./data/DNN_miRNA.pth",
        help="SiCmiR model")
    args = parser.parse_args()

    predict(
        model_path=args.model,
        input_path=args.input,
        output_path=args.output,
        output_dir=args.output_dir,
        storage_dir=args.storage_dir,
        extract=args.extract,
        normalization=args.normalization,
        save_extract=args.save_extract,       
        save_zscore_input=args.save_zscore_input
    )

if __name__ == "__main__":
    main()

