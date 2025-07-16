# SiCmiR
**SiCmiR** predicts miRNA expression profile from only 977 LINCS L1000 landmark genes!

---
## Installation
You can install and run **SiCmiR** using either `pip` or `conda`.

**Option 1:  Install with pip**
```bash
pip install -r Requirements.txt
```
**Option 2: Create Conda environment (CPU version)**
```bash
conda create -n SiCmiR python=3.8 pandas numpy scipy pytorch -c pytorch -y
```
###Option 3: Use GPU (Optional)
If your system supports CUDA and you want to use GPU acceleration, install PyTorch with the appropriate CUDA version. Visit:
👉 https://pytorch.org/get-started/locally/

```bash
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
```
**Then activate your environment:**
```bash
conda activate SiCmiR
```
