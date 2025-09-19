import scanpy as sc
import matplotlib.pyplot as plt


def perform_dimensionality_reduction(adata):
    """Perform PCA, neighbor computation, and UMAP for clustering."""
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.9, flavor="igraph", n_iterations=2, directed=False)
    return adata


def identify_marker_genes(adata):
    """Identify and visualize marker genes for each cluster."""
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, save="_ranked_genes.png")
    plt.close()

    # Extract top genes per cluster
    cluster_top = {
        group: adata.uns["rank_genes_groups"]["names"][group][:10]
        for group in adata.uns["rank_genes_groups"]["names"].dtype.names
    }
    return cluster_top


def annotate_clusters(adata):
    """Annotate clusters with cell type labels."""
    cluster_annotations = {
        "0": "T cells (memory/activated)",
        "1": "T cells (naïve/effector)",
        "2": "B cells",
        "3": "Monocytes / DC",
        "4": "NK cells",
        "5": "Monocytes (classical)",
        "6": "Antigen-presenting cells (DC-like)",
        "7": "Megakaryocytes",
        "8": "Proliferating cells (cycling)",
    }
    adata.obs["cell_type"] = (
        adata.obs["leiden"].map(cluster_annotations).astype("category")
    )
    return adata, cluster_annotations
