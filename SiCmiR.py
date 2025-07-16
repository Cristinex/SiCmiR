# coding=utf-8
import sys
import argparse
#----Part I: Data preparation----#
import pandas as pd
import numpy as np
from scipy.stats import zscore


parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", default="L1000_mRNA.csv")
parser.add_argument("--output", "-o", default="result.csv")
parser.add_argument("--storage_dir", "-sd", default='./data')
parser.add_argument("--extract", "-e", default=False)
parser.add_argument("--normalization", "-norm", default=False)
parser.add_argument("--RNA", "-r", default='extracted_zscored_mRNA.csv')
parser.add_argument("--extracted_output", "-eo",default = 'extracted_mRNA.csv')
parser.add_argument("--output_dir", "-od", default='./')
args = parser.parse_args()

def build_submatrix(df_all: pd.DataFrame, picked_names, picked_alias, fill_value=0, zscore_norm = False):
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
    if zscore_norm:
        sub = sub.apply(lambda x: zscore(x, ddof=0), axis=0).fillna(0) #zscore scaling and fill NA with 0
        print('normalized')
    return sub

#----Part II: Prediction----#
import torch
import torch.nn as nn

from torch.nn import init
from torch.utils.data import DataLoader
torch.set_default_dtype(torch.float)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

class NeuralNet(nn.Module):
    def __init__(self, input_size=977, hidden_size=1024, output_size=1298):
        super(NeuralNet, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.ReLU(),
            nn.BatchNorm1d(hidden_size),
            nn.Dropout(p=0.3),  # dropout训
            nn.Linear(hidden_size, hidden_size),
            
            #nn.ReLU(),
            #nn.BatchNorm1d(hidden_size),
            #nn.Dropout(p=0.3),
            
            nn.ReLU(),
            nn.BatchNorm1d(hidden_size),
            nn.Dropout(p=0.3),
            
            nn.Linear(hidden_size, output_size),
        )
    def forward(self, x):
        x = self.model(x)
        return x

# define a loss function
loss_fn = torch.nn.MSELoss().to(device)
net=torch.load(f"{args.storage_dir}/DNN_miRNA.pth",map_location=device)
net.eval()


if args.extract:
    L1000_gene = pd.read_csv(f'{args.storage_dir}/L1000gene.csv',header=None) #Genes required as input
    picked_names = L1000_gene.loc[:,0]    # L1000 genes in pre-defined nomenclature
    picked_alias = L1000_gene.loc[:,1]    # gene names alias 
    data_matrix = pd.read_csv(args.input,index_col=0).T #!!!input matrix here with genes in rows and samples in columns!!!

    extracted_zscored_mRNA = build_submatrix(data_matrix,picked_names, picked_alias,zscore_norm=args.normalization)  
    
    extracted_zscored_mRNA.to_csv(args.extracted_output)

    X_test_P = extracted_zscored_mRNA 
else:
    X_test_P = pd.read_csv(args.input,index_col=0)
X_test_Tensor = torch.tensor(X_test_P.values,dtype=torch.float)
miRNA_1298 = pd.read_csv(f"{args.storage_dir}/1298miRNA.csv", header = None).loc[:,0]
y_test_out = pd.DataFrame(net(X_test_Tensor).detach().numpy())
y_test_out.columns = miRNA_1298
y_test_out.index = X_test_P.index

print(y_test_out)
y_test_out.to_csv(f'{args.output_dir}/{args.output}') 





