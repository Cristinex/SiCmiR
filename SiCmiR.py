# coding=utf-8
import os
import argparse
import pandas as pd
import numpy as np
from scipy.stats import zscore
import torch
import torch.nn as nn
import logging
import random
import math
import time
import scanpy as sc

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

torch.set_default_dtype(torch.float)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

#----Part I: Data preparation----#
def build_submatrix(df_all: pd.DataFrame, picked_names, picked_alias, output_dir, fill_value=0, save_extract=None):
    logger.info("Building submatrix with 977 landmark genes.")
    sub = pd.DataFrame({
    name: df_all.get(name, df_all.get(name_alias, pd.Series(fill_value, index=df_all.index)))
    for name, name_alias in zip(picked_names, picked_alias)})
    if save_extract:
        
        sub.to_csv(os.path.join(output_dir, save_extract))
        logger.info(f"Extracted unzscored submatrix saved to {os.path.join(output_dir, save_extract)}")
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

    logger.info("Starting SiCmiR prediction pipeline.")
    if not os.path.exists(input_path):
        logger.error(f"Input file {input_path} does not exist.")
        raise FileNotFoundError(f"Input file {input_path} not found.")
    if not os.path.isdir(storage_dir):
        logger.error(f"Input file {storage_dir} does not exist.")
        raise FileNotFoundError(f"Input file {storage_dir} not found.")
    L1000_file = os.path.join(storage_dir, 'L1000gene.csv')
    miRNA_file = os.path.join(storage_dir, '1298miRNA.csv')
    model_file = os.path.join(storage_dir, 'DNN_miRNA.pth')
    if not os.path.exists(L1000_file):
        logger.error(f"L1000 gene file {L1000_file} does not exist.")
        raise FileNotFoundError(f"L1000 gene file {L1000_file} not found.")
    if not os.path.exists(miRNA_file):
        logger.error(f"miRNA file {miRNA_file} does not exist.")
        raise FileNotFoundError(f"miRNA file {miRNA_file} not found.")
    if not os.path.exists(model_file):
        logger.error(f"miRNA file {model_file} does not exist.")
        raise FileNotFoundError(f"miRNA file {model_file} not found.")        
    if not os.path.isdir(output_dir) and '.':
        logger.warning(f"Input file {output_dir} does not exist. Created.")
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    if len(output_path) == 0:
        logger.error(f"Expected output filename.")
        raise ValueError("No output filename specified.")


    if input_path.endswith(".h5ad"):
        ad = sc.read_h5ad(input_path)
        data_matrix = ad.to_df()
    elif input_path.endswith((".csv", ".csv.gz")):
        data_matrix = pd.read_csv(input_path, index_col=0).T
    else: 
        logger.error(f"Invalid input formmat.")
        raise TypeError ("Unsupported file format; expected .h5ad or .csv(.gz)")
    logger.info(f"Input matrix shape: {data_matrix.shape}")
    if data_matrix.empty:
        logger.error("Input matrix is empty or improperly formatted.")
        raise ValueError("Invalid input matrix.")
        

    if extract:
        L1000_gene = pd.read_csv(os.path.join(storage_dir,'L1000gene.csv'), header=None)
        picked_names = L1000_gene.loc[:, 0]
        picked_alias = L1000_gene.loc[:, 1]
        extracted_mRNA = build_submatrix(
            data_matrix, picked_names, picked_alias, output_dir, save_extract=save_extract)
    else:
        extracted_mRNA = data_matrix

    if normalization:
        extracted_zscored_mRNA = extracted_mRNA.apply(lambda x: zscore(x, ddof=0), axis=0).fillna(0) #zscore scaling and fill NA with 0
        if save_zscore_input:
            extracted_zscored_mRNA.to_csv(os.path.join(output_dir,save_zscore_input))
        X_test_P = extracted_zscored_mRNA
    else:
        X_test_P = extracted_mRNA

    # Load model
    net = torch.load(f"{model_path}", map_location=device)
    net.eval()
    
    X_test_Tensor = torch.tensor(X_test_P.values, dtype=torch.float)
    miRNA_1298 = pd.read_csv(os.path.join(storage_dir,'1298miRNA.csv'), header=None).loc[:, 0]
    y_test_out = pd.DataFrame(net(X_test_Tensor).detach().numpy())
    y_test_out.columns = miRNA_1298
    y_test_out.index = X_test_P.index

    # Save predicted miRNA expression
    y_test_out.to_csv(os.path.join(output_dir,output_path))
    logger.info(f"Results saved to {os.path.join(output_dir, output_path)}")
    
    return y_test_out

def main():
    parser = argparse.ArgumentParser(description="SiCmiR: Predict miRNA expression from mRNA data")
    parser.add_argument("--input", "-i", type=str, default="./demo/test_mRNA.csv",
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

