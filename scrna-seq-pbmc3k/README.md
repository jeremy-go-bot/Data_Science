🧬 scRNA-seq Analysis: PBMC3k Dataset
🚀 Overview
This project showcases an advanced single-cell RNA sequencing (scRNA-seq) pipeline applied to the PBMC3k dataset (3,000 peripheral blood mononuclear cells).
By integrating bioinformatics and machine learning, I demonstrate expertise in data preprocessing, clustering, marker gene identification, visualization, and predictive modeling to uncover biological insights.
Goal: Build a robust, reproducible pipeline to identify and classify immune cell types with precision.

🔬 Key Features
* Data Processing: Loaded and preprocessed 10x Genomics data using Scanpy, with filtering, normalization, and highly variable gene selection.
* Clustering: Applied PCA, UMAP, and Leiden clustering to identify 9 distinct immune cell clusters (e.g., T cells, B cells, NK cells).
* Marker Genes: Identified top marker genes (e.g., CD79A, NKG7) using Wilcoxon rank-sum test for accurate cell type annotation.
* Machine Learning: Trained a Random Forest classifier (AUC: 0.98) to predict cell types, with SHAP analysis to explain feature importance.
* Visualizations: Generated professional UMAP plots, marker gene heatmaps, and SHAP summary plots.

📊 Results
* Clusters: Identified and annotated 9 immune cell types with high confidence.
* Model Performance: Achieved a Random Forest AUC of 0.98, showcasing robust classification.
* Key Genes: Highlighted critical markers (e.g., CD79A for B cells, NKG7 for NK cells) driving biological insights.
* Visualizations: Produced clear, publication-quality plots to visualize clusters and gene expression.

🛠️ Setup and Usage
1. Clone the Repository

git clone https://github.com/jeremy-go-bot/Data_Science.git
cd Data_Science/scrna-seq-pbmc3k
2. Install Dependencies

pip install -r requirements.txt
Required packages: scanpy, matplotlib, seaborn, scikit-learn, shap, pandas, numpy.
3. Download Data
Get the PBMC3k dataset from 10x Genomics Extract into:

data/filtered_gene_bc_matrices/hg19/
4. Run the Pipeline

python scrna_seq_analysis.py
5. Outputs (saved in figures/)
* umap_clusters.png, umap_predicted_cell_types.png
* marker_heatmap.png
* shap_summary_plot.png
* cluster_annotations_and_top_genes.md

📁 Artifacts
* Code: scrna_seq_analysis.ipynb (modular, documented pipeline)
* Figures: UMAP, heatmap, and SHAP plots in figures/
* Report: Cluster annotations and top genes in figures/cluster_annotations_and_top_genes.md

🎯 Why This Matters
This project demonstrates my ability to:
* Build end-to-end bioinformatics pipelines with Scanpy.
* Apply machine learning for precise cell type classification.
* Create publication-quality visualizations and interpretable models.
* Deliver reproducible, well-documented code for real-world impact.
Explore the code and report to see the pipeline in action!

📚 Data Source
10x Genomics – PBMC3k Dataset
