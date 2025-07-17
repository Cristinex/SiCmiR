# SiCmiR

**SiCmiR** predicts miRNA expression profile from only 977 LINCS L1000 landmark genes!

Visit **SiCmiR-Atlas** for database: https://awi.cuhk.edu.cn/~SiCmiR/

---

## Usage
### 💡 Script: 
**SiCmiR.py** support extract landmark genes, conduct z-score normalization & miRNA prediction (the same as script/SiCmiR_beta.py)

**SiCmiR_full.py** support generate pseudo bulk samples, extract landmark genes, conduct z-score normalization & miRNA prediction

The main script accepts input gene expression matrix (genes in rows and samples/cells in columns) and outputs predicted miRNA expression.

### 📂 Directory tree
```bash
SiCmiR-main/
├── SiCmiR.py                   # main script
├── Requirements.txt            # dependency version
├── README.md                   # README

├── data/                       # supporting documents for main script
│   ├── L1000gene.csv           # landmark genes (used in extraction)
│   ├── miRNA1298.csv           # miRNAs (used for miRNA prediction)
│   └── DNN_miRNA.pth           # pretrained_model (used for miRNA prediction)

├── script/                     # script
│   ├── SiCmiR_beta.py          # the same as SiCmiR.py
│   └── SiCmiR_full.py          # ⚡ copy it to project_root (SiCmiR-main) before use

├── demo/                       # data for quick test (genes in rows and samples/cells in columns)
│   ├── test_mRNA.csv           # an mRNA expression matrix already extracted & normalized and no need to generate pseudo bulk samples
│   │── raw_GSE64465_mRNA.h5ad   # a single-cell expression matrix from Seurat after batch effect removal containing 3k+ cells
│   │                             (See download link at raw_GSE64465_mRNA.txt)
│   └── cell_meta.csv           # cell metadata for GSE644565

└── outputs/                    # demonstration for output
│   └── predicted_miRNA.csv     # demonstration for output
│   └── predicted_GSE64465_miRNA.csv     # demonstration for output
```

### 📝 Example

### Example 1: Predict using a preprocessed matrix 
(already extracted & normalized and no need to generate pseudo bulk samples)
```bash
python SiCmiR.py --input ./demo/test_mRNA.csv --output predicted_miRNA.csv
```
### Example 2: Use full pipeline 
(generate pseudo bulk samples, extract landmark genes, conduct z-score normalization & miRNA prediction)
```bash
python SiCmiR_full.py \
  --input ./demo/raw_GSE64465_mRNA.h5ad \
  --output_dir ./ \
  --save_extract extracted_unzscore_mRNA.csv \
  --save_zscore_input extracted_zscore_mRNA.csv \
  --output predicted_GSE64465_miRNA.csv\
  --normalization\
  --pooling_method bootstrap \
  --group_file ./demo/cell_meta.csv
```

### 📕 Available Arguments
| Argument           | Alias       | Default                          | Description |
|--------------------|-------------|----------------------------------|-------------|
| `--input`          | `-i`        | `'./demo/test_mRNA.csv'`      | Input mRNA expression matrix <br>.csv: genes in rows, samples in columns; <br> .h5ad: genes in columns, samples in rows;<br>📌 Remove batch effects before use if needed. |
| `--output`         | `-o`        | `'predicted_miRNA.csv'`          | Output file: predicted miRNA expression, with miRNAs in rows and samples in columns |
| `--output_dir`     | `-od`       | `'./'`                            | Output directory for saving the predicted results |
| `--storage_dir`    | `-sd`       | `'./data/'`                       | Directory to store pretrained model and auxiliary files |
| `--no_extract`     | `-n`        | `False`                           | Do not extract 977 landmark genes from input mRNA matrix (Not Recomended) |
| `--normalization`  | `-norm `    | `False`                           | Whether to perform z-score normalization on input matrix.<br>⚠️ Avoid repeating normalization |
| `--save_extract`   | `-se`       | `None`                            | Output file name for extracted but unnormalized mRNA matrix (if `--extract` is used) |
| `--save_zscore_input` | `-z`     | `None`                            | Save zscored mRNA matrix used directly for prediction (if `--normalization` is used) |
| `--model`          | `-m`        | `'./data/DNN_miRNA.pth'`          | Path to SiCmiR model |
| `--pooling_method`          | `-p`        | `None`          | Choose a pooling method. choices=[None, 'average', 'bootstrap']  |
| `--group_file`          | `-g`        | `'./demo/group.csv'`          | CSV file with sample groups (index: sample/cell names, column 'group': group labels) |
| `--random_seed`          | `-r`        | `system entropy`          | Random seed for bootstrap sampling. If not provided, uses system entropy (Not recommend to set) |
---
## Requiremnt Installation &  Download

### 📌 Setup the environment using either `pip` or `conda`.

**Option 1:**  Install with pip
```bash
pip install -r Requirements.txt
```
**Option 2:** Create Conda environment (CPU version)
```bash
conda create -n SiCmiR python=3.8 pandas numpy scipy pytorch -c pytorch -y

or

conda create -n SiCmiR python=3.8 pandas=1.5 numpy=1.21 scipy=1.10 pytorch=2.0 -c pytorch -y
```

**Then activate your environment:**
```bash
conda activate SiCmiR
```
**Option 3:** Use GPU (Optional & Not Recommended since there is no training process applying SiCmiR)

If your system supports CUDA and you want to use GPU acceleration, install PyTorch with the appropriate CUDA version. Visit:
https://pytorch.org/get-started/locally/

```bash
conda create -n SiCmiR python=3.8 pandas numpy scipy
conda activate SiCmiR
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
```
### 📌 Obtain SiCmiR in two ways:

**Option 1: Clone the repository**

If you have `git` installed and internet access to GitHub:

```bash
git clone https://github.com/Cristine/SiCmiR.git
cd SiCmiR
python SiCmiR.py --input ./demo/test_mRNA.csv --output predicted_miRNA.csv
```
This method allows you to keep the project up to date via git pull.

**Option 2: Download as ZIP**

If you don't use git:

Go to https://github.com/Cristine/SiCmiR

- Click the green Code button → Download ZIP
- Unzip the downloaded SiCmiR-main.zip & Navigate to the extracted folder:
```bash
unzip SiCmiR-main.zip
cd SiCmiR-main
python SiCmiR.py --input ./demo/test_mRNA.csv --output predicted_miRNA.csv
```




