import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import os
from src.preprocess import load_and_preprocess_data
from src.analysis import (
    perform_dimensionality_reduction,
    identify_marker_genes,
    annotate_clusters,
)
from src.visualization import visualize_clusters
from src.modeling import train_random_forest, explain_model

st.title("🧬 scRNA-seq Analysis Dashboard: PBMC3k")
st.markdown(
    "Explore single-cell RNA sequencing analysis results for the PBMC3k dataset."
)

# Load and preprocess data
data_path = st.text_input(
    "Enter path to PBMC3k dataset", "data/filtered_gene_bc_matrices/hg19"
)
if st.button("Run Analysis"):
    with st.spinner("Processing data..."):
        adata = load_and_preprocess_data(data_path)
        adata = perform_dimensionality_reduction(adata)
        cluster_top = identify_marker_genes(adata)
        adata, cluster_annotations = annotate_clusters(adata)

        # Display UMAP plot
        st.subheader("UMAP Visualization")
        sc.pl.umap(adata, color=["cell_type"], save="_clusters.png")
        st.image(os.path.join("figures", "umap_clusters.png"))

        # Display marker gene table
        st.subheader("Cluster Annotations and Top Marker Genes")
        df = pd.DataFrame(
            [
                [cluster, cluster_annotations.get(cluster, "Unknown"), ", ".join(genes)]
                for cluster, genes in cluster_top.items()
            ],
            columns=["Cluster", "Cell Type", "Top Genes"],
        )
        st.table(df)

        # Display SHAP plot
        st.subheader("Random Forest Feature Importance (SHAP)")
        markers = {
            "T cells (memory/activated)": [
                "LTB",
                "IL32",
                "HINT1",
                "TRAF3IP3",
                "SELL",
                "GIMAP5",
                "CD2",
                "GIMAP7",
                "LDLRAP1",
                "MAL",
            ],
            "T cells (naïve/effector)": [
                "LTB",
                "IL32",
                "CD2",
                "AQP3",
                "GIMAP7",
                "HINT1",
                "GIMAP5",
                "MAL",
                "FYB",
                "GIMAP4",
            ],
            "B cells": [
                "CD79A",
                "CD79B",
                "HLA-DPB1",
                "MS4A1",
                "HLA-DQA1",
                "HLA-DRB1",
                "CD37",
                "HLA-DQB1",
                "HLA-DPA1",
                "TCL1A",
            ],
            "Monocytes / DC": [
                "LST1",
                "COTL1",
                "FCER1G",
                "AIF1",
                "IFITM3",
                "SAT1",
                "CTSS",
                "S100A11",
                "FCGR3A",
                "CST3",
            ],
            "NK cells": [
                "NKG7",
                "CST7",
                "GZMA",
                "CTSW",
                "CCL5",
                "PRF1",
                "GZMM",
                "PTPRCAP",
                "GZMB",
                "FGFBP2",
            ],
            "Monocytes (classical)": [
                "S100A8",
                "TYROBP",
                "FCN1",
                "S100A6",
                "CST3",
                "LGALS2",
                "LGALS1",
                "GPX1",
                "GSTP1",
                "AIF1",
            ],
            "Antigen-presenting cells (DC-like)": [
                "HLA-DPA1",
                "HLA-DPB1",
                "HLA-DRB1",
                "HLA-DQA1",
                "CST3",
                "FCER1A",
                "HLA-DMA",
                "HLA-DQB1",
                "CLEC10A",
                "GRN",
            ],
            "Megakaryocytes": [
                "PPBP",
                "GPX1",
                "SDPR",
                "PF4",
                "NRGN",
                "CCL5",
                "GNG11",
                "SPARC",
                "RGS18",
                "HIST1H2AC",
            ],
            "Proliferating cells (cycling)": [
                "KIAA0101",
                "PCNA",
                "CKS1B",
                "BIRC5",
                "STMN1",
                "HPRT1",
                "SNRPB",
                "PSMD14",
                "STOML2",
                "MCM3",
            ],
        }
        palette = {
            "T cells (memory/activated)": "#1f77b4",
            "T cells (naïve/effector)": "#ff7f0e",
            "B cells": "#2ca02c",
            "Monocytes / DC": "#d62728",
            "NK cells": "#9467bd",
            "Monocytes (classical)": "#8c564b",
            "Antigen-presenting cells (DC-like)": "#e377c2",
            "Megakaryocytes": "#7f7f7f",
            "Proliferating cells (cycling)": "#bcbd22",
        }
        rf, X_test, marker_genes = train_random_forest(adata, markers, palette)
        explain_model(rf, X_test, marker_genes)
        st.image(os.path.join("figures", "shap_summary_plot.png"))

st.markdown(
    "Source dataset: [10x Genomics](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)"
)
st.markdown("Source code: [GitHub](https://github.com/jeremy-go-bot/Data_Science)")
