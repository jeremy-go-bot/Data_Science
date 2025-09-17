```
# 🧬 scRNA-seq Analysis: PBMC3k Dataset
# Single-cell RNA sequencing pipeline for PBMC3k dataset (3,000 PBMCs).
# Integrates bioinformatics and ML for preprocessing, clustering, marker genes, and predictions.

## 🚀 Overview
# Goal: Build a robust pipeline to classify immune cell types with precision.
# Showcases preprocessing, clustering, marker gene analysis, and predictive modeling.

## 🔬 Key Features
# - Data: Preprocess 10x Genomics data with Scanpy (filter, normalize, HVG selection).
# - Clustering: PCA, UMAP, Leiden to identify 9 immune cell clusters (T cells, B cells, NK cells).
# - Markers: Wilcoxon rank-sum test to find markers (e.g., CD79A, NKG7).
# - ML: Random Forest classifier (AUC: 0.98) with SHAP for feature importance.
# - Visuals: UMAP plots, marker heatmaps, SHAP summary plots.

## 📊 Results
# - Clusters: 9 immune cell types annotated with high confidence.
# - Model: Random Forest AUC = 0.98.
# - Genes: CD79A (B cells), NKG7 (NK cells) as key markers.
# - Visuals: Publication-quality plots.

## 🛠️ Setup
# 1. Clone Repo
#    git clone https://github.com/jeremy-go-bot/Data_Science.git
#    cd Data_Science/scrna-seq-pbmc3k
# 2. Install Dependencies
#    pip install -r requirements.txt
#    # Needs: scanpy, matplotlib, seaborn, scikit-learn, shap, pandas, numpy
# 3. Download Data
#    # Get PBMC3k from 10x Genomics: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
#    # Extract to: data/filtered_gene_bc_matrices/hg19/
# 4. Run Pipeline
#    python scrna_seq_analysis.py
# 5. Outputs (in figures/)
#    - umap_clusters.png, umap_predicted_cell_types.png
#    - marker_heatmap.png, shap_summary_plot.png
#    - cluster_annotations_and_top_genes.md

## 📁 Artifacts
# - Code: scrna_seq_analysis.ipynb
# - Figures: UMAP, heatmaps, SHAP in figures/
# - Report: Annotations and genes in figures/cluster_annotations_and_top_genes.md

## 📚 Data Source
# 10x Genomics – PBMC3k Dataset
```