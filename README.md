# SiCmiR
**SiCmiR** predicts miRNA expression profile from only 977 LINCS L1000 landmark genes!

---

## Usage
Script: SiCmiR.py

The main script accepts input gene expression data (977 landmark genes) and outputs predicted miRNA expression.

## Example

### Example 1: Predict using a preprocessed matrix (already extracted & normalized)
```bash
python SiCmiR.py --input ./demo/test_mRNA.csv --output predicted_miRNA.csv
```
### Example 2: Use default pretrained model + full pipeline
```bash
python SiCmiR.py \
  --input ./demo/raw_test_mRNA.csv \
  --normalization\
  --save_extract extracted_unzscore_mRNA.csv \
  --save_zscore_input 'extracted_zscore_mRNA.csv' \
  --output_dir ./ \
  --output predicted_miRNA.csv
```

### Available Arguments
| Argument           | Alias       | Default                          | Description |
|--------------------|-------------|----------------------------------|-------------|
| `--input`          | `-i`        | `'./example/test_mRNA.csv'`      | Input mRNA expression matrix (genes in rows, samples in columns).<br>📌 Remove batch effects before use if needed. |
| `--output`         | `-o`        | `'predicted_miRNA.csv'`          | Output file: predicted miRNA expression, with miRNAs in rows and samples in columns. |
| `--output_dir`     | `-od`       | `'./'`                            | Output directory for saving the predicted results. |
| `--storage_dir`    | `-sd`       | `'./data/'`                       | Directory to store pretrained model and auxiliary files. |
| `--no_extract`     | `-n`        | `False`                           | Do not extract 977 landmark genes from input mRNA matrix (Not Recomended). |
| `--normalization`  | `-norm `    | `False`                           | Whether to perform z-score normalization on input matrix.<br>⚠️ Avoid repeating normalization. |
| `--save_extract`   | `-se`       | `None`                            | Output file name for extracted but unnormalized mRNA matrix (if `--extract` is used) |
| `--save_zscore_input` | `-z`     | `None`                            | Save zscored mRNA matrix used directly for prediction (if `--normalization` is used). |
| `--model`          | `-m`        | `'./data/DNN_miRNA.pth'`          | Path to SiCmiR model. |


---
## Requiremnt Installation &  Download

You can install and run **SiCmiR** using either `pip` or `conda`.

**Option 1:**  Install with pip
```bash
pip install -r Requirements.txt
```
**Option 2:** Create Conda environment (CPU version)
```bash
conda create -n SiCmiR python=3.8 pandas numpy scipy pytorch -c pytorch -y
```
**Option 3:** Use GPU (Optional & Not Recommended since there is no training process in applying SiCmiR)

If your system supports CUDA and you want to use GPU acceleration, install PyTorch with the appropriate CUDA version. Visit:
https://pytorch.org/get-started/locally/

```bash
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
```
**Then activate your environment:**
```bash
conda activate SiCmiR
```

You can obtain SiCmiR in two ways:

### Option 1: Clone the repository (recommended)

If you have `git` installed and internet access to GitHub:

```bash
git clone https://github.com/Cristine/SiCmiR.git
cd SiCmiR
python SiCmiR.py --input ./example/test_mRNA.csv --output predicted_miRNA.csv
```
This method allows you to keep the project up to date via git pull.

### Option 2: Download as ZIP

If you don't use git:

Go to https://github.com/Cristine/SiCmiR

- Click the green Code button → Download ZIP
- Unzip the downloaded SiCmiR-main.zip
- Navigate to the extracted folder:
```bash
cd SiCmiR-main
python SiCmiR.py --input ./example/test_mRNA.csv --output predicted_miRNA.csv
```




