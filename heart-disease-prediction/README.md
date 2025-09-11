# ❤️ Heart Disease Prediction - UCI Dataset

A complete **classification modeling project** using the [UCI Heart Disease (Cleveland Clinic) dataset](https://www.kaggle.com/datasets/cherngs/heart-disease-cleveland-uci).  
This project walks through the full ML workflow: **EDA, preprocessing, feature engineering, model building, hyperparameter tuning, and explainability with SHAP**.

---

## 🎯 Objective
The aim is to **predict the presence of heart disease** (binary target: `condition`, 1 = disease, 0 = no disease) based on medical and demographic patient attributes.  

Goals:
- Identify the most influential risk factors.  
- Build and evaluate several machine learning models.  
- Balance accuracy with **interpretability for clinical use**.  

---

## 📊 Workflow

1. **Data Inspection**
   - Dataset: **297 rows × 14 columns**.  
   - Features:  
     - Numerical → `age`, `trestbps`, `chol`, `thalach`, `oldpeak`  
     - Categorical → `sex`, `cp`, `fbs`, `restecg`, `exang`, `slope`, `ca`, `thal`  
   - Class balance: 46% positive cases, 54% negative.  
   - No missing values or duplicates.  
   - Male-dominant dataset: 67% men, 33% women.  

2. **Exploratory Data Analysis (EDA)**
   - Heart disease patients showed:  
     - Lower max heart rate (`thalach`).  
     - Higher ST depression (`oldpeak`).  
     - Stronger chest pain symptoms (`cp`).  
   - Key predictors:  
     - Positive correlation → `thal`, `ca`, `cp`.  
     - Negative correlation → `thalach`.  
   - Multicollinearity among `oldpeak`, `slope`, `exang`.  

3. **Preprocessing & Feature Engineering**
   - Created engineered features:  
     - `oldpeak_log` = log-transformed ST depression.  
     - `trestbps_age` = resting blood pressure × age.  
     - `stress_score` = `oldpeak_log × slope` (synthetic stress indicator).  
   - Scaling: `RobustScaler` for numerical features.  
   - One-hot encoding for categorical variables.  
   - Final dataset: **297 × 20 features**.  

4. **Modeling**
   - **Baseline**: Logistic Regression (ROC AUC ≈ 0.914) outperformed tree-based models.  
   - **After Feature Engineering**:  
     - Logistic Regression ROC AUC ≈ 0.918  
     - SVC ROC AUC ≈ 0.900  
   - **Final Model (Optuna-tuned Logistic Regression):**  
     - Train ROC AUC = 0.919  
     - Test ROC AUC = 0.975  
     - Average Precision = 0.979  
     - F1 Score = 0.902  

---

## 🏆 Results

| Model                | ROC AUC | Avg Precision | F1   |
|----------------------|---------|---------------|------|
| Logistic Regression  | 0.918   | 0.918         | 0.820 |
| Random Forest        | 0.885   | 0.883         | 0.796 |
| XGBoost              | 0.870   | 0.871         | 0.773 |
| LightGBM             | 0.885   | 0.882         | 0.771 |
| SVC                  | 0.900   | 0.895         | 0.813 |

👉 **Best Model: Logistic Regression (Optuna-tuned)** — excellent performance and high interpretability, making it suitable for clinical contexts.  

---

## 🔍 Key Insights
- Strongest predictors of heart disease:  
  - **Chest pain type (`cp`)** ↑  
  - **Number of major vessels (`ca`)** ↑  
  - **Thalassemia test result (`thal`)** ↑  
  - **Stress score (oldpeak × slope)** ↑  
  - **Maximum heart rate (`thalach`)** ↓  
- **Key lesson**: On small medical datasets, **feature engineering + interpretable models** often outperform complex models while preserving trust in clinical settings.  

---

## 📂 Repository Structure

📦 heart-disease-prediction
┣ 📜 heart_disease_notebook.ipynb # Full workflow
┣ 📜 heart_disease_model.pkl # Final Logistic Regression pipeline
┣ 📂 figures # Plots and visualizations
┣ 📜 README.md # Project documentation


---

## 🚀 Future Work
- Explore **stacking ensembles** (LogReg + SVC + tree-based).  
- Perform **cross-cohort validation** with other UCI subsets.  
- Conduct **fairness analysis** (e.g., gender imbalance).  

---


