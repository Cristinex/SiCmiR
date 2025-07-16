# SiCmiR
**SiCmiR** predicts miRNA expression profile from only 977 LINCS L1000 landmark genes!

---
## Installation
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

##Usage
Script: SiCmiR.py
The main script accepts input gene expression data (977 landmark genes) and outputs predicted miRNA expression.

Example
```bash
python SiCmiR.py --input ./example/example_input.csv --output predicted_miRNA.csv

Available arguments
Argument	Alias	Default	Description
--input	-i	./example/test_mRNA.csv	Input mRNA expression matrix (genes in rows, samples in columns).
⚠️ Remove batch effects before use if needed.
--output	-o	./example/result.csv	Output file: predicted miRNA expression, with miRNAs in rows and samples in columns.
--output_dir	-od	'./'	Output directory for saving the predicted results.
--storage_dir	-sd	'./data/'	Directory to store pretrained model and auxiliary files.
--extract	-e	False	Whether to extract 977 landmark genes from a full mRNA matrix.
--normalization	-norm	False	Whether to perform z-score normalization on input matrix.
⚠️ Avoid repeating normalization.
--RNA	-r	'extracted_zscored_mRNA.csv'	Input RNA matrix after extraction & normalization, used directly for prediction.
--extracted_output	-eo	'extracted_unzscored_mRNA.csv'	Output file for extracted but unnormalized mRNA matrix (if --extract True).

