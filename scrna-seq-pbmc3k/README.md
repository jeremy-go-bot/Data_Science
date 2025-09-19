# 🧫 scRNA-seq Analysis: PBMC3k Dataset

## 🚀 Overview
This project showcases an advanced **single-cell RNA sequencing (scRNA-seq)** pipeline applied to the **PBMC3k dataset** (3,000 peripheral blood mononuclear cells).  

By integrating **bioinformatics**, **machine learning**, and a **Streamlit dashboard**, I demonstrate expertise in:  
- Data preprocessing  
- Clustering  
- Marker gene identification  
- Visualization  
- Predictive modeling  
- Reproducibility via **Docker containerization**  

**Goal:** Build a robust, reproducible pipeline to identify and classify immune cell types with precision.

---

## 🔬 Key Features
- **Data Processing**: Loaded and preprocessed 10x Genomics data using Scanpy (filtering, normalization, highly variable gene selection).  
- **Clustering**: Applied PCA, UMAP, Leiden clustering → identified 9 distinct immune cell clusters (e.g., T cells, B cells, NK cells).  
- **Marker Genes**: Identified top markers (e.g., CD79A, NKG7) using Wilcoxon rank-sum test.  
- **Machine Learning**: Trained a Random Forest classifier (AUC = 0.98) with SHAP analysis for feature importance.  
- **Visualizations**: Produced UMAP plots, heatmaps, SHAP plots (publication-quality).  
- **Dashboard**: Interactive Streamlit app displaying UMAP clusters, marker tables, SHAP feature importance.  
- **Reproducibility**: Fully containerized with Docker for easy deployment and portability. 

---

## 📊 Results

| **Aspect**            | **Outcome** |
|------------------------|-------------|
| **Clusters**          | 9 immune cell types annotated with high confidence |
| **Model Performance** | Random Forest → AUC = 0.98 |
| **Key Genes**         | CD79A (B cells), NKG7 (NK cells) |
| **Visualizations**    | UMAP plots, marker heatmaps, SHAP summaries |

---

## 🖥️ Streamlit Dashboard

Run the interactive dashboard and explore:  
- **UMAP Visualization** of clusters  
- **Cluster annotations & marker genes**  
- **SHAP feature importance** of Random Forest model  

---

## 🛠️ Setup and Usage

### Option 1 – Local (pip/conda)

### 1. Clone the Repository
git clone https://github.com/jeremy-go-bot/Data_Science.git  
cd Data_Science/scrna-seq-pbmc3k  

### 2. Install Dependencies
pip install -r requirements.txt  
or
conda env create -f environment.yml

Required packages: scanpy, matplotlib, seaborn, scikit-learn, shap, pandas, numpy , streamlit

### 3. Download Data
Download the PBMC3k dataset from 10x Genomics (https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).

Extract into:  
data/filtered_gene_bc_matrices/hg19/  

### 4. Run the Pipeline
streamlit run app.py

### 5. Outputs (saved in figures/)
- umap_clusters.png  
- umap_predicted_cell_types.png  
- marker_heatmap.png  
- shap_summary_plot.png  
- cluster_annotations_and_top_genes.md  

### Option 2 – Run with Docker (recommended 🚀)
Build and run the containerized app:

docker build -t scrna-seq-app

docker run -p 8501:8501 scrna-seq-app

Then access the dashboard at: http://localhost:8501

💡 If port 8501 is already in use, you can map to another port:


docker run -p 8502:8501 scrna-seq-app
Then open http://localhost:8502

---

## 📁 Artifacts
- **Code**: scrna_seq_analysis.ipynb (modular, documented pipeline, src/) + Streamlit app (app.py)  
- **Figures**: UMAP, heatmap, SHAP plots in figures/  
- **Report**: Cluster annotations & top genes in figures/cluster_annotations_and_top_genes.md  
- **Dockerfile**: Ensures full reproducibility and portability
---

## 🎯 Why This Matters
This project demonstrates my ability to:
- Build end-to-end bioinformatics pipelines with Scanpy  
- Apply machine learning for precise cell type classification  
- Create publication-quality visualizations and interpretable models  
- Deliver reproducible, well-documented code for real-world impact  
- Ensure reproducibility with Docker (one-command setup, no pip/conda conflicts)


---

## 📚 Data Source
10x Genomics – PBMC3k Dataset (https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).
