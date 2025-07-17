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
def pool_samples(df: pd.DataFrame, group_col=None, pooling_method='none', bootstrap_n=2, bootstrap_fraction=0.8, random_seed=42):
    if pooling_method == 'none' or group_col is None:
        logger.info("No pooling applied.")
        return df
    
    if not isinstance(group_col, pd.Series):
        if len(group_col) != len(df.index):
            logger.error(f"Length of group_col ({len(group_col)}) does not match number of rows in input. {len(df.index)} expected.")
            raise ValueError("Group column length mismatch.")
        group_col = pd.Series(group_col, index=df.index)
    
    unique_groups = group_col.unique()
    pooled_dfs = []
    
    if pooling_method == 'average':
        logger.info("Pooling samples by averaging.")
        for group in unique_groups:
            group_samples = group_col[group_col == group].index
            if len(group_samples) > 0:
                pooled_df = df.loc[group_samples].mean(axis=0,numeric_only=True).to_frame(name=group).T
                pooled_dfs.append(pooled_df)
    elif pooling_method == 'bootstrap':
        random.seed(random_seed)
        for group in unique_groups:
            group_samples = group_col[group_col == group].index
            bootstrap_n = max(1, math.ceil(len(group_samples) / 100))
            logger.info(f"Pooling {group} samples/cells by bootstrap sampling with {bootstrap_n} iterations and {bootstrap_fraction*100}% non-replacing sampling. This may take time. Wait patiently and do not stop the code, please.")

            if len(group_samples) > 0:
                for i in range(bootstrap_n):
                    sampled_cols = np.random.choice(group_samples, size=int(len(group_samples) * bootstrap_fraction), replace=False)
                    pooled_df = df.loc[sampled_cols].mean(axis=0,numeric_only=True).to_frame(name=f"{group}_bootstrap_{i+1}").T
                    pooled_dfs.append(pooled_df)
    
    if pooled_dfs:
        pooled_df = pd.concat(pooled_dfs, axis=0)
        logger.info(f"Pooled DataFrame shape: {pooled_df.shape}")
        if pooled_df.shape[0] <= 2:
            logger.warning("Invalid number of pseudo_bulk sample. Returning original DataFrame")
            return df
        return pooled_df
    else:
        logger.warning("No valid groups found for pooling. Returning original DataFrame.")
        return df

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
            nn.Dropout(p=0.3),
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
    save_zscore_input=None,
    pooling_method=None,
    group_file=None,
    bootstrap_n=2,
    random_seed=None
):

    logger.info("Starting SiCmiR prediction pipeline.")
    if random_seed is None:
        random_seed = int(time.time())
        logger.info(f"Using dynamic seed: {random_seed}")
    else:
        logger.info(f"Using provided seed: {random_seed}")
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
    if group_file and not os.path.exists(group_file):
        logger.error(f"Group file {group_file} does not exist.")
        raise FileNotFoundError(f"Group file {group_file} not found.")
    if len(output_path) == 0:
        logger.error(f"Expected output filename.")
        raise ValueError("No output filename specified.")


    # Load and preprocess input data

    
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
    # Load group information for pooling (if provided)
    group_col = None
    if group_file and pooling_method:
        if pooling_method not in ['average','bootstrap']:
            logger.warning('Invalid pooling method. Skip pooling.')
            pooling_method = None
        else:
            group_df = pd.read_csv(group_file, index_col=0)
            if 'group' in group_df.columns:
                group_col = group_df['group']
                if group_col is not None:
                    group_col = group_col.reindex(data_matrix.index)
                    missing_samples = group_col.index[group_col.isna()]
                    missing_in_group = [s for s in data_matrix.index if s not in group_col.index]
                if missing_samples.tolist():
                    logger.warning(f"Missing group labels for samples: {missing_samples.tolist()}. These samples will be ignored in pooling.")
                extra_samples = [s for s in group_col.index if s not in data_matrix.index]
                if extra_samples:
                    logger.warning(f"Group file contains extra samples not in input matrix: {extra_samples}. These will be ignored.")
                    group_col = group_col[group_col.index.isin(data_matrix.index)]
                logger.info(f"Loaded group information from {group_file}")
                if missing_in_group:
                    logger.warning(f"Samples in input matrix not found in group file: {missing_in_group}. These will be excluded from pooling.")
                    data_matrix = data_matrix.loc[data_matrix.index.isin(group_col.index)]
            else:
                logger.error("Group file does not contain a column named 'group'. Required for pooling.")
                raise ValueError("Invalid group file: 'group' column missing.")
    
    # Pool samples
    pooled_matrix = pool_samples(
        data_matrix,
        group_col=group_col,
        pooling_method=pooling_method,
        bootstrap_n=bootstrap_n,
        bootstrap_fraction=0.8,
        random_seed=random_seed
    )
    logger.info(f"Pooled matrix shape: {pooled_matrix.shape}")


    # Extract 977 landmark genes
    if extract:
        L1000_gene = pd.read_csv(os.path.join(storage_dir, 'L1000gene.csv'), header=None)
        picked_names = L1000_gene.loc[:, 0]
        picked_alias = L1000_gene.loc[:, 1]
        extracted_mRNA = build_submatrix(
            pooled_matrix, picked_names, picked_alias, output_dir, save_extract=save_extract)
    else:
        extracted_mRNA = pooled_matrix
    
    # Normalize if requested
    if normalization:
        logger.info("Applying z-score normalization.")
        extracted_zscored_mRNA = extracted_mRNA.apply(lambda x: zscore(x, ddof=0), axis=0).fillna(0)
        if save_zscore_input:
            extracted_zscored_mRNA.to_csv(os.path.join(output_dir, save_zscore_input))
            logger.info(f"Z-scored matrix saved to {os.path.join(output_dir, save_zscore_input)}")
        X_test_P = extracted_zscored_mRNA
    else:
        X_test_P = extracted_mRNA
    
    # Load model and predict
    logger.info(f"Loading model from {model_path}")
    net = torch.load(model_path, map_location=device)
    net.eval()
    
    X_test_Tensor = torch.tensor(X_test_P.values, dtype=torch.float)
    miRNA_1298 = pd.read_csv(os.path.join(storage_dir, '1298miRNA.csv'), header=None).loc[:, 0]
    y_test_out = pd.DataFrame(net(X_test_Tensor).detach().numpy())
    y_test_out.columns = miRNA_1298
    y_test_out.index = X_test_P.index
    
    # Save predicted miRNA expression
    y_test_out.to_csv(os.path.join(output_dir, output_path))
    logger.info(f"Results saved to {os.path.join(output_dir, output_path)}")
    
    return y_test_out

def main():
    parser = argparse.ArgumentParser(description="SiCmiR: Predict miRNA expression from mRNA data")
    parser.add_argument("--input", "-i", type=str, default="./demo/test_mRNA.csv",
        help="Input mRNA matrix (genes in rows, samples in columns). Default: %(default)s")
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
        help="Perform z-score normalization on input matrix. Default: False")
    parser.add_argument("--save_extract", "-se", default=None,
        help="Output file name for extracted but unnormalized mRNA matrix. Default: None")
    parser.add_argument("--save_zscore_input", "-z", default=None,
        help="Save z-scored mRNA matrix used for prediction. Default: None")
    parser.add_argument("--model", "-m", default="./data/DNN_miRNA.pth",
        help="SiCmiR model")
    parser.add_argument("--pooling_method", "-p", choices=[None, 'average', 'bootstrap'], default='none',
        help="Pooling method: 'none', 'average', or 'bootstrap'. Default: %(default)s")
    parser.add_argument("--group_file", "-g", default=None,
        help="CSV file with sample groups (index: sample names, column 'group': group labels). Default: None")
    parser.add_argument("--random_seed", "-r", type=int, default=None,
        help="Random seed for bootstrap sampling. If not provided, uses system entropy. Default: None")
    
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
        save_zscore_input=args.save_zscore_input,
        pooling_method=args.pooling_method,
        group_file=args.group_file,
        random_seed=args.random_seed
    )

if __name__ == "__main__":
    main()
