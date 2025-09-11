# 🦠 COVID-19 Open Research Dataset Challenge (CORD-19 NLP)

An end-to-end **Natural Language Processing project** using the [CORD-19 dataset](https://www.kaggle.com/datasets/allen-institute-for-ai/CORD-19-research-challenge) (Allen Institute / Kaggle).  
This project demonstrates how to: **clean large-scale corpora, engineer features, perform topic modeling, and build document retrieval systems with both classical and transformer-based NLP methods**.

---

## 🎯 Objective
The goal is to analyze and organize COVID-19 scientific literature to:  
- Identify **temporal trends** and thematic shifts in publications.  
- Extract research themes via **topic modeling (LDA)**.  
- Build a **document retrieval system** using TF-IDF and cosine similarity.  
- Leverage **BioBERT embeddings** for semantic search beyond keyword matching.  

---

## 📊 Workflow

1. **Data Loading & Inspection**
   - ~800,000 papers with metadata: title, abstract, authors, journal, publish_time.  
   - Full-text parses from JSON (`pdf_json`).  
   - Missing values >50% for some fields (authors, journal, abstracts).  

2. **Data Cleaning**
   - Removed duplicates (`cord_uid`).  
   - Dropped missing titles, imputed `publish_time` with median.  
   - Normalized author/journal fields.  
   - Filtered abstracts to 50–5000 characters.  
   - COVID-era filter (≥ Dec 2019).  
   - Final dataset: **~600,000 clean records**.  

3. **Exploratory Data Analysis (EDA)**
   - Publication spike after 2020.  
   - Top sources: WHO, Medline, PMC.  
   - Top journals: *PLoS One*, *bioRxiv*, *medRxiv*.  
   - >200k unique authors → global collaboration.  
   - Word cloud highlighted frequent terms: *COVID-19, pandemic, SARS-CoV-2*.  

4. **Feature Engineering**
   - Temporal features: year, month, COVID-era flag.  
   - Keyword flags: `has_risk_factor`, `has_transmission`, `has_vaccine`.  
   - Abstract representation with **TF-IDF (5,000 features)**.  
   - Full-text extracted with **Dask** for scalability.  
   - Stored in **Parquet** format for efficiency.  

5. **Topic Modeling (LDA)**
   - 10 latent topics extracted from abstracts.  
   - Example topics:  
     - *Epidemiology*: covid, cases, transmission, population  
     - *Clinical*: patients, treatment, hospital, symptoms  
     - *Vaccines*: vaccine, efficacy, immune, trial  
   - Perplexity & coherence validated interpretability.  

6. **Document Retrieval**
   - TF-IDF + cosine similarity enabled **keyword-based search**.  
   - Example: *“What are known risk factors for COVID-19?”* → studies on comorbidities, mortality, demographics.  

7. **Semantic Search with BioBERT**
   - Pretrained model: `dmis-lab/biobert-v1.1`.  
   - Dense embeddings (768D) generated via mean pooling.  
   - Cosine similarity for retrieval.  
   - Outperformed TF-IDF by retrieving semantically related results (e.g., *“risk factors”* → hypertension, diabetes, obesity).  

---

## 🏆 Key Learnings
- Large-scale biomedical corpora require **rigorous cleaning & normalization**.  
- **TF-IDF** is strong for fast keyword retrieval but limited by vocabulary.  
- **LDA** provides interpretable thematic clusters for bibliometrics.  
- **BioBERT embeddings** enable semantic search beyond keyword matching.  
- Hybrid systems (**TF-IDF + BioBERT**) balance **speed and accuracy** in retrieval.  

---

## 📂 Repository Structure
📦 covid19-cord19-nlp
┣ 📜 Covid19_notebook.ipynb # Full workflow
┣ 📂 models # Saved models
┣ 📂 figures # Visualizations and plots
┣ 📜 README.md # Project documentation

---
