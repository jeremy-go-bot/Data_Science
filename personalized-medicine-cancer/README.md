# 🧬 Personalized Medicine: Redefining Cancer Treatment

A Kaggle competition project ([MSKCC](https://www.kaggle.com/c/msk-redefining-cancer-treatment)) demonstrating **hybrid tabular + NLP modeling** for cancer mutation classification.  
This project integrates structured genomic features with unstructured biomedical text, combining **XGBoost + transformer embeddings** in a unified pipeline.

---

## 🎯 Objective
The goal is to **predict the clinical impact of cancer mutations** (multi-class classification, 9 classes) by leveraging both **tabular genomic features (Gene, Variation)** and **long clinical abstracts**.  

Key highlights:
- Handle imbalanced multi-class classification.  
- Combine tabular + NLP features in one pipeline.  
- Use transformer embeddings (**MiniLM**) for efficient, high-quality biomedical text representation.  
- Ensure interpretability with **SHAP** to validate feature importance and biological insights.  
- Build a **deployment-ready pipeline** with Kaggle submission output.  

---

## 📊 Workflow

1. **Data Loading & Cleaning**
   - Merged `training_variants` (tabular) with `training_text` (abstracts).  
   - Cleaned gene/variation strings, standardized case.  
   - Text preprocessing: lowercasing, punctuation removal, stopword filtering.  
   - Integrity checks: no duplicate IDs, missing abstracts handled.  

2. **Exploratory Data Analysis (EDA)**
   - Class imbalance: Class 7 ≈ 30% of samples.  
   - Most frequent genes: **TP53, BRCA1**.  
   - Mutation types: **missense mutations dominate**; deletions/insertions rarer.  
   - Abstract length: ~3k–5k tokens.  
   - Gene–class dependencies observed (heatmaps).  

   📂 Figures: `class_distribution.png`, `top_genes.png`, `gene_class_heatmap.png`, `text_length_tokens.png`, `variation_type_by_class.png`  

3. **Feature Engineering**
   - **Tabular Features**:  
     - `Gene_Freq` (frequency encoding)  
     - `Variation_Type` (one-hot)  
     - `Variation_Pos` (numeric extraction)  
     - `Var_Freq` (frequency encoding)  
   - **Text Features**:  
     - TF-IDF (5k bi-grams) → baseline  
     - **MiniLM embeddings (384-dim)** → transformer encoder, efficient and accurate  
   - Final set: Tabular + MiniLM (~400 features).  

4. **Modeling & Evaluation**
   - Pipeline: **SMOTE → XGBoost**.  
   - Cross-validation (5-fold Stratified):  
     - CV Log Loss ≈ 1.10  
     - CV F1 Macro ≈ 0.51  
   - Hyperparameter tuning with **Optuna**.  
   - Explainability with **SHAP**:  
     - MiniLM embeddings dominate.  
     - Tabular features add complementary signal.  

   📂 Figures: `shap_importance.png`, `confusion_matrix.png`, `precision_recall_curves.png`, `ablation_study.png`  

5. **Ensemble Strategy**
   - TF-IDF baseline model built.  
   - Averaged predictions with MiniLM pipeline → improved robustness.  
   - Ensemble Log Loss: ~1.21 (vs. 1.10 for MiniLM).  

6. **Inference & Submission**
   - Trained final pipeline on full dataset.  
   - Generated `submission.csv` with required Kaggle format.  
   - Deployment-ready workflow.  

---

## 🏆 Results

| Model          | CV Log Loss |
|----------------|-------------|
| TF-IDF + XGB   | ~1.39       |
| MiniLM + XGB   | ~1.10       |
| Ensemble       | ~1.21       |

| Features       | CV Log Loss | Std Dev |
|----------------|-------------|---------|
| Full           | 1.10        | 0.034   |
| No MiniLM      | 1.44        | 0.040   |
| No Tabular     | 1.20        | 0.041   |
| Only Gene Freq | 1.62        | 0.041   |

👉 **Best Model: MiniLM + XGB (Optuna tuned)**  
- Achieved the strongest balance of accuracy and interpretability.  
- Outperformed classical baselines (e.g., Logistic Regression, TF-IDF).  

---

## 🔍 Key Learnings
- **Hybrid features** (tabular + text) outperform single-modality models.  
- **MiniLM embeddings** provide powerful yet efficient text representations.  
- **SHAP** validates domain insights: TP53, BRCA1, and missense mutations are key.  
- Systematic experiments (CV, Optuna, ablation) ensure **robust and reproducible results**.  

---

## 📂 Repository Structure
📦 personalized-medicine-cancer
┣ 📜 personalized_medicine_notebook.ipynb # Full workflow
┣ 📂 figures # Visualizations and plots
┣ 📜 README.md # Project documentation

---
