import scanpy as sc


def load_and_preprocess_data(data_path):
    """Load and preprocess scRNA-seq data from 10x Genomics format."""
    adata = sc.read_10x_mtx(data_path, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()

    # Filter low-quality cells and genes
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Select highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]

    # Scale data
    sc.pp.scale(adata, max_value=10)

    return adata
