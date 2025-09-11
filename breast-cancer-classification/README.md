# 🧬 Breast Cancer Classification (Wisconsin Diagnostic)

A complete **end-to-end Machine Learning project** using the [Breast Cancer Wisconsin dataset](https://www.kaggle.com/datasets/uciml/breast-cancer-wisconsin-data) (UCI / scikit-learn).  
This project covers the full ML pipeline: **EDA, preprocessing, feature analysis, dimensionality reduction, model building, hyperparameter tuning, and evaluation**.

---

## 🎯 Objective
The goal is to **predict whether a tumor is malignant or benign** based on 30 numerical features extracted from digitized images of breast cell nuclei.  

Main aims:
- Identify the tumor characteristics most associated with malignancy.  
- Compare linear and non-linear classifiers on biomedical tabular data.  
- Build a reproducible, high-performing pipeline with automated hyperparameter optimization.  

---

## 📊 Workflow

1. **Data Inspection**
   - 569 samples, 33 columns (30 features, 1 target, 2 irrelevant).  
   - Dropped `id` and `Unnamed: 32`.  
   - No missing or duplicate values.  

2. **Exploratory Data Analysis (EDA)**
   - Target distribution: 212 malignant (37%), 357 benign (63%) → imbalanced.  
   - Malignant tumors had higher values for radius, perimeter, area, concave points.  
   - Correlation analysis showed strong multicollinearity among size-related features.  
   - Variance Inflation Factor (VIF) confirmed redundancy → dimensionality reduction required.  

3. **Preprocessing**
   - **PowerTransformer + StandardScaler**: normalized skewed features.  
   - **PCA (n_components=0.95):** reduced 30 → ~10 components, preserving 95% variance.  
   - **SMOTE** applied (only on training folds) to balance classes.  

4. **Modeling & Hyperparameter Tuning**
   - Models: Logistic Regression, Random Forest, XGBoost, LightGBM, SVC.  
   - Hyperparameters optimized with **Optuna**.  
   - Experiment tracking via **MLflow**.  

5. **Evaluation**
   - Metrics: ROC AUC, Average Precision (AP), F1 Score.  
   - Stratified K-Fold CV for balanced splits.  

---

## 🏆 Results

| Model                | ROC AUC | Avg Precision | F1 Score |
|-----------------------|---------|---------------|----------|
| Logistic Regression   | 0.995   | 0.994         | 0.964    |
| Random Forest         | 0.987   | 0.983         | 0.932    |
| XGBoost               | 0.992   | 0.989         | 0.944    |
| LightGBM              | 0.990   | 0.986         | 0.923    |
| SVC                   | 0.994   | 0.993         | 0.961    |

👉 **Best Model: LightGBM (Optuna-tuned)**  
- Hyperparameters: `n_estimators=157`, `learning_rate≈0.29`, `max_depth=3`  
- Final Test Results: **ROC AUC = 0.997 | AP = 0.996 | F1 = 0.963**  

---

## 🔍 Key Learnings
- **PCA** reduced dimensionality and mitigated multicollinearity without sacrificing variance.  
- **SMOTE** balanced the dataset, preventing bias toward benign cases.  
- **Tree-based methods** (LightGBM, XGBoost) captured nonlinear feature interactions better than linear models.  
- A systematic pipeline with preprocessing + automated tuning delivers **robust, reproducible results**.  

---

## 📂 Repository Structure
📦 breast-cancer-classification
┣ 📜 breast_cancer_notebook.ipynb # Full workflow
┣ 📂 figures # Plots and visualizations
┣ 📜 README.md # Project documentation

---
