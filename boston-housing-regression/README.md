# 🏠 Boston Housing Prices - Advanced Regression Techniques

A complete **regression modeling project** using the [Boston Housing dataset](https://www.kaggle.com/datasets/fedesoriano/the-boston-houseprice-data).  
This project demonstrates the full machine learning pipeline: **EDA, preprocessing, feature engineering, baseline regression models, advanced tree-based methods, and neural networks**.

---

## 🎯 Objective
The aim is to **predict housing prices (MEDV)** — the median value of owner-occupied homes (in $1000s) — based on socio-economic, demographic, and environmental features.  

Goals:
- Explore factors influencing housing prices.  
- Compare different regression models, from linear to nonlinear.  
- Achieve the best trade-off between accuracy, interpretability, and generalization.  

---

## 📊 Workflow

1. **Data Inspection**
   - Dataset: **506 rows × 14 columns**  
   - Features: numeric, with `CHAS` binary and `RAD` categorical-like.  
   - No missing values; `CRIM` highly skewed.  

2. **Exploratory Data Analysis (EDA)**
   - **Target (MEDV)** capped at 50 (artificial ceiling).  
   - Strong correlations:  
     - `RM` (avg. rooms) ↑ → higher MEDV.  
     - `LSTAT` (% lower-status) ↑ → lower MEDV.  
   - Wealth indicators (`TAX`, `PTRATIO`, `NOX`, `INDUS`) clustered together.  
   - Outliers in `CRIM` and capped MEDV identified.  

3. **Preprocessing**
   - Log transforms applied to `CRIM`, `ZN`, `LSTAT`.  
   - `RobustScaler` used (resistant to outliers).  
   - Dropped feature `B` (ethical issues, weak signal).  
   - Built preprocessing pipeline with `ColumnTransformer`.  

4. **Feature Engineering**
   - **Mutual Information** confirmed `RM` and `LSTAT` as strongest predictors.  
   - Added **polynomial features** (quadratic + interactions).  
   - Tested **PCA** for dimensionality reduction after expansion.  

5. **Modeling**
   - **Linear Regression (baseline):**  
     - RMSE ≈ 3.90 | MAE ≈ 2.48 | R² ≈ 0.79  
   - **XGBoost (Optuna-tuned):**  
     - RMSE ≈ 2.31 | MAE ≈ 1.72 | R² ≈ 0.93  
     - Best performer, capturing nonlinear effects.  
     - SHAP confirmed importance of `LSTAT`, `RM`, `AGE`, `DIS`.  
   - **Neural Network (PyTorch MLP):**  
     - RMSE ≈ 4.87 | R² ≈ 0.67  
     - Underperformed due to small dataset. Early stopping applied.  

---

## 🏆 Results

| Model                 | RMSE | MAE  | R²   |
|------------------------|------|------|------|
| Linear Regression      | 3.90 | 2.48 | 0.79 |
| XGBoost (Optuna-tuned) | 2.31 | 1.72 | 0.93 |
| Neural Network (MLP)   | 4.87 | 3.41 | 0.67 |

👉 **Best Model: XGBoost**, combining high accuracy with interpretability through SHAP analysis.  

---

## 🔍 Key Insights
- Housing prices increase with **larger homes (RM ↑)** and decrease with **higher % of lower-status population (LSTAT ↑)**.  
- Wealth-related features cluster together, influencing MEDV.  
- **Tree-based methods outperform linear and neural networks** on structured tabular data.  
- Neural networks are not always effective on small datasets.  

---

## 📂 Repository Structure

📦 boston-housing-regression
┣ 📜 boston_housing_notebook.ipynb # Full workflow
┣ 📜 boston_house_price_model.joblib # Final XGBoost model + preprocessing
┣ 📂 figures # Generated figures
┣ 📜 README.md # Project documentation


---

## 🚀 Future Work
- Add **regularized linear models** (Ridge, Lasso, ElasticNet).  
- Explore **stacking/blending ensembles**.  
- Consider removing capped MEDV values to improve generalization.  

---
