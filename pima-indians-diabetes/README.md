# 🩺 Pima Indians Diabetes Prediction

A complete **end-to-end Machine Learning project** using the [Pima Indians Diabetes dataset](https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database) (UCI / Kaggle).  
This project demonstrates the full workflow: **EDA, feature engineering, preprocessing, model building, hyperparameter tuning, and interpretability with SHAP**.

—

## 🎯 Objective
The goal is to **predict whether a patient has diabetes** (binary classification: 0 = no diabetes, 1 = diabetes) using 8 medical attributes.  

Main aims:
- Identify the clinical factors most associated with diabetes.  
- Compare linear and non-linear models on biomedical data.  
- Build a reproducible, high-performing ML pipeline with automated tuning.  
- Provide **model interpretability (SHAP)** to align predictions with clinical knowledge.  

---

## 📊 Workflow

1. **Data Inspection**
   - Dataset: 768 samples, 9 columns (8 features + 1 target).  
   - No missing or duplicate values.  
   - Target imbalance: ~34% positive cases (diabetes).  

2. **Exploratory Data Analysis (EDA)**
   - Invalid medical values: many zeros in *Glucose, BloodPressure, SkinThickness, Insulin, BMI*.  
   - Outliers in *Insulin* and *DiabetesPedigreeFunction*.  
   - Correlations:  
     - Strongest predictors → `Glucose` (0.47), `BMI` (0.31), `Age` (0.24).  
     - `Pregnancies` correlated with `Age` (0.61).  

3. **Feature Engineering**
   - Imputed invalid zeros using **median by class**.  
   - Applied **RobustScaler** to reduce outlier influence.  
   - Skewness analysis flagged *Insulin* & *DiabetesPedigreeFunction*.  

4. **Modeling & Hyperparameter Tuning**
   - Models tested: Logistic Regression, Random Forest, XGBoost, LightGBM, CatBoost, SVC.  
   - Class imbalance addressed with **SMOTE** (applied only on training folds).  
   - Hyperparameter optimization with **Optuna** (RF & CatBoost).  
   - Ensemble: **Stacking (RF + CatBoost)** tested for boost.  

5. **Evaluation**
   - Metrics: ROC AUC, Average Precision (AP), F1, Recall.  
   - Stratified K-Fold CV ensured balanced evaluation.  

---

## 🏆 Results

| Model              | ROC AUC | Avg Precision | F1 Score |
|--------------------|---------|---------------|----------|
| Logistic Regression| 0.833   | 0.824         | 0.713    |
| Random Forest      | 0.889   | 0.879         | 0.811    |
| XGBoost            | 0.864   | 0.845         | 0.789    |
| LightGBM           | 0.872   | 0.844         | 0.792    |
| CatBoost           | **0.885** | **0.864**   | **0.802** |
| SVC                | 0.869   | 0.849         | 0.790    |
| Stacking (RF+CB)   | **0.901** | **0.888**   | **0.828** |

👉 **Best Model (single): CatBoost (Optuna-tuned)**  
- Hyperparameters: `depth = 10`, `l2_leaf_reg ≈ 0.0019`  
- Final Test Performance: **ROC AUC = 0.912 | AP = 0.888 | F1 = 0.824**  

👉 **Best Overall: Stacking (RF + CatBoost)**  
- Slightly higher metrics, but less interpretable → **CatBoost chosen as final model** for clinical relevance.  

---

## 🔍 Key Learnings
- **Zero imputation + robust scaling** improved data quality and performance.  
- **Tree-based models** (RF, CatBoost, LightGBM) outperformed linear baselines.  
- **Optuna hyperparameter tuning** pushed CatBoost to near-optimal results.  
- **SHAP interpretation** confirmed medical validity:  
  - High **Glucose** & **BMI** → higher diabetes risk.  
  - **Age** also a strong predictor.  
- Interpretability is crucial for clinical adoption.  

---

## 📂 Repository Structure
📦 pima-indians-diabetes
┣ 📜 pima_indian_diabetes_notebook.ipynb # Full workflow
┣ 📂 figures # Plots and visualizations
┣ 📜 README.md # Project documentation

---
