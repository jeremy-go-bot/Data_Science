import scanpy as sc
import matplotlib.pyplot as plt


def visualize_clusters(adata, palette):
    """Visualize UMAP with annotated cell types."""
    sc.pl.umap(adata, color=["cell_type"], palette=palette, save="_clusters.png")
    plt.close()
