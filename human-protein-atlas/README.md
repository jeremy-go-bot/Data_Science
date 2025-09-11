# 🧬 Human Protein Atlas Image Classification (HPAIC)

A deep learning project based on the [Human Protein Atlas Image Classification challenge](https://www.kaggle.com/c/human-protein-atlas-image-classification).  
This project demonstrates how to adapt **Convolutional Neural Networks (CNNs)** for **multi-label protein classification** on biomedical fluorescence microscopy images.  

---

## 🎯 Objective
The task is to **classify proteins in human cells** from microscopy images:  
- Each sample includes **4 channels** (Red, Green, Blue, Yellow), representing different protein stains.  
- Images may belong to **multiple categories simultaneously** (28 protein localization classes).  

---

## 📊 Workflow

1. **Data Loading & Preprocessing**
   - Dataset: `train.csv` + ~30k microscopy images.  
   - Labels converted to **multi-label one-hot encoding** across 28 classes.  
   - Severe class imbalance detected (rare categories underrepresented).  
   - Built a custom **PyTorch Dataset** to handle 4-channel inputs.  

2. **Baseline CNN**
   - Implemented a simple 3-layer CNN.  
   - Input shape: `(4, 224, 224)` tensors.  
   - Output: 28 sigmoid-activated probabilities.  
   - Loss: **Binary Cross-Entropy (BCE)**.  
   - Achieved convergence, but underfit complex patterns.  

3. **Transfer Learning (ResNet18)**
   - Adapted **ResNet18 (pretrained on ImageNet)**:  
     - Modified first convolution to accept 4 channels.  
     - Replaced fully connected layer with 28-class output.  
   - Training with **BCEWithLogitsLoss** stabilized learning.  
   - Metrics: **Macro F1, Micro F1**.  
   - Outperformed the baseline CNN significantly.  

---

## 🏆 Results

| Model        | Macro F1 | Micro F1 |
|--------------|----------|----------|
| Simple CNN   | ~0.052   | ~0.347   |
| ResNet18 (TL)| ~0.212   | ~0.504   |

👉 Even with minimal tuning, **ResNet18 adapted to 4 channels** showed clear improvements over the baseline CNN.  

---

## 🔍 Key Learnings
- **Multi-label tasks** require sigmoid outputs + BCE loss (instead of softmax).  
- Pretrained CNNs can be **easily adapted** to handle non-standard channel inputs.  
- Class imbalance is a major challenge → strategies like **Focal Loss, re-weighting, or oversampling** are needed.  
- Stronger architectures (EfficientNet, ConvNeXt) and advanced augmentation could push performance further.  

---

## 📂 Repository Structure
📦 human-protein-atlas
┣ 📜 human_protein_notebook.ipynb # Full workflow
┣ 📂 figures # Visualizations and plots
┣ 📜 README.md # Project documentation

---
