# 🧬 Gene Expression Cancer RNA-Seq Classification

A complete **end-to-end Machine Learning project** using the [Gene Expression Cancer RNA-Seq dataset](https://www.kaggle.com/datasets/manjitbaishya001/gene-expression-cancer-rna-seq).  
This project covers: **EDA, dimensionality reduction, feature selection, modeling, hyperparameter tuning, and gene-level interpretation with SHAP**.

---

## 🎯 Objective
The goal is to **classify cancer types** from RNA-Seq gene expression profiles (~20,000 features).  

Specific aims:
- Identify transcriptomic signatures distinguishing cancer subtypes.  
- Apply **dimensionality reduction (PCA, t-SNE)** to explore high-dimensional structure.  
- Compare linear and nonlinear classifiers on large-scale biomedical data.  
- Derive **biologically meaningful gene signatures** with SHAP values.  

---

## 📊 Workflow

1. **Data Inspection**
   - 801 samples, 20,531 gene features + 1 target.  
   - No missing or duplicate values.  
   - Memory optimized (float64 → float32).  

2. **Exploratory Data Analysis (EDA)**
   - Mild class imbalance (entropy ≈ 0.94).  
   - Expression distributions: many zeros → log1p considered.  
   - Skewness: mean ≈ 0.92 (positively skewed).  
   - Strong multicollinearity among genes → motivated PCA.  

3. **Dimensionality Reduction**
   - **PCA (2D):** ~19% variance explained, partial subtype separation.  
   - **t-SNE (2D):** clear clustering, preserved local structure.  
   - **PCA (80% variance):** 168 components retained.  

4. **Preprocessing**
   - Standard scaling.  
   - Variance thresholding to remove uninformative genes.  
   - PCA applied → best generalization performance.  

5. **Modeling & Evaluation**
   - Models tested: Logistic Regression, Random Forest, XGBoost, LightGBM, CatBoost, SVC.  
   - Metrics: **F1 Macro, Recall Macro, Balanced Accuracy** with Stratified CV.  
   - Results (with PCA, 80% variance explained):  

| Model               | F1 Macro | Recall Macro | Balanced Acc |
|---------------------|----------|--------------|--------------|
| Logistic Regression | **1.000**| **1.000**    | **1.000**    |
| Random Forest       | 0.982    | 0.977        | 0.977        |
| XGBoost             | 0.970    | 0.968        | 0.968        |
| SVC                 | 0.991    | 0.988        | 0.988        |
| LightGBM            | 0.980    | 0.976        | 0.976        |
| CatBoost            | 0.978    | 0.975        | 0.975        |

👉 **Best Model: Logistic Regression + PCA**  
- Achieved **perfect classification (F1 = 1.0)**.  
- Simpler, interpretable, and computationally efficient compared to tree ensembles.  

6. **Gene-Level Interpretation**
   - Combined **VarianceThreshold + Random Forest + SHAP** to identify top 100 predictive genes.  
   - Each gene contributed ~0.5–1.5% to classification.  
   - Multi-gene signatures aligned with biological expectations (cancer is polygenic).  
   - Exported list: `top_100_genes_shap.csv` → usable for pathway enrichment analysis.  

---

## 🏆 Key Learnings
- High-dimensional RNA-Seq requires **dimensionality reduction** to handle noise and collinearity.  
- **PCA + Logistic Regression** surprisingly outperformed tree-based ensembles.  
- **t-SNE** provided intuitive clustering of cancer subtypes in transcriptomic space.  
- SHAP enabled biologically relevant interpretation, yielding a **gene signature**.  
- Robust preprocessing allows **simple models** to match or outperform complex ones.  

---

## 📂 Repository Structure
📦 gene-expression-cancer
┣ 📜 RNA-SEQ_notebook.ipynb # Full workflow
┣ 📜 top_100_genes_shap.csv # Extracted top 100 predictive genes
┣ 📂 figures # Visualizations and plots
┣ 📜 README.md # Project documentation

---
