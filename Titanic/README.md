# 🚢 Titanic Survival Prediction

A complete **end-to-end Machine Learning project** using the [Titanic dataset](https://www.kaggle.com/datasets/yasserh/titanic-dataset) from Kaggle.  
This project demonstrates the full data science workflow: **data exploration, preprocessing, feature engineering, modeling, hyperparameter tuning, and ensembling**.

---

## 🎯 Objective
The goal is to **predict passenger survival on the Titanic** based on demographic and travel information (e.g., age, sex, ticket class, fare, family size).  

Key aims:
- Identify the most influential survival factors.  
- Build and compare multiple machine learning models.  
- Optimize accuracy, interpretability, and generalization for Kaggle competition settings.  

---

## 📊 Workflow

1. **Data Exploration & Cleaning**
   - Inspected datasets (Train: 891 rows, Test: 418 rows).  
   - Handled missing values (`Age`, `Embarked`, `Fare`), dropped `Cabin`.  
   - Visualized survival rates by gender, class, and fare distribution.

2. **Feature Engineering**
   - Created new features:
     - `FamilySize` (siblings/spouses + parents/children)  
     - `Title` (extracted from passenger names)  
     - `FareLog` (log-transformed fare to reduce skew)  
   - One-hot encoded categorical variables.  

3. **Modeling & Tuning**
   - **Baseline**: Logistic Regression (~84% accuracy).  
   - **Advanced models**:  
     - XGBoost (tuned with Optuna)  
     - LightGBM (tuned with Optuna)  
     - Stacking Ensemble (LogReg + XGB + LGBM).  

4. **Evaluation**
   - Compared models using Accuracy, Precision, Recall, F1, and ROC AUC.  
   - Plotted ROC curves for performance visualization.  

---

## 🏆 Results

| Model                | Accuracy | Precision | Recall | F1 Score | ROC AUC |
|-----------------------|----------|-----------|--------|----------|---------|
| Logistic Regression   | 0.84     | 0.82      | 0.75   | 0.79     | 0.83    |
| XGBoost               | 0.83     | 0.80      | 0.74   | 0.77     | 0.84    |
| LightGBM              | 0.83     | 0.79      | 0.77   | 0.78     | 0.83    |
| **Stacking Ensemble** | **0.85** | **0.84**  | **0.77** | **0.80** | **0.87** |

👉 **Best Model: Stacking Ensemble**, offering the most balanced performance.  

---

## 🔍 Key Insights
- Women, higher-class passengers, and those paying higher fares had the highest survival probability.  
- Feature engineering (`FamilySize`, `Title`, `FareLog`) improved model performance significantly.  
- Ensemble methods combined the strengths of linear and nonlinear models for better generalization.  

---

## 📂 Repository Structure

📦 titanic-survival-prediction
┣ 📜 titanic_notebook.ipynb # Full workflow notebook
┣ 📜 figures # All figures
┣ 📜 README.md # Project documentation


---

## 🚀 Future Work
- Test alternative ensemble strategies (e.g., blending, soft voting).  
- Apply SMOTE to balance survival vs. non-survival classes.  
- Explore deep learning models for comparison.  

---
