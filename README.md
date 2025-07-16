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
--input", "-i", default="./example/test_mRNA.csv": input mRNA matrix with genes in rows and samples in columns
                                                  (Removal of batch effects before use is preferable).
--output_dir", "-od", default='./': output directory for predicted miRNA expression profile
--output", "-o", default="./example/result.csv": output predicted miRNA expression profiles with miRNA in
                                                 rows and samples in columns.
--storage_dir", "-sd", default='./data/': directory for storage of model parameters and supportive files.
--extract", "-e", default=False: whether need to execute extraction from original mRNA expression matrix.
--normalization", "-norm", default=False: whether need to conduct zscore normalization for the
                                          input/extracted mRNA expression matrix (should avoid repeatly zscore)
--RNA", "-r", default='extracted_zscored_mRNA.csv': input mRNA expression profile for model
--extracted_output", "-eo",default = 'extracted_unzscored_mRNA.csv': unzscored mRNA expression
                                                                      profiles extracted if -e True
