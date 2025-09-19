from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import shap
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os

os.environ["FIGURES_DIR"] = "./figures"
FIGURES_DIR = os.environ["FIGURES_DIR"]


def train_random_forest(adata, markers, palette):
    """Train and evaluate a Random Forest classifier for cell type prediction."""
    marker_genes = [
        gene for genes in markers.values() for gene in genes if gene in adata.var_names
    ]
    X = adata[:, marker_genes].X
    y = adata.obs["cell_type"].values

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(X_train, y_train)

    # Predict and evaluate
    y_prob = rf.predict_proba(X_test)
    auc = roc_auc_score(y_test, y_prob, multi_class="ovr")
    print(f"Random Forest AUC: {auc:.4f}")

    # Predict on full dataset for visualization
    adata.obs["predicted_cell_type"] = rf.predict(adata[:, marker_genes].X)
    sc.pl.umap(
        adata,
        color=["predicted_cell_type"],
        palette=palette,
        save="_predicted_cell_types.png",
    )
    plt.close()

    return rf, X_test, marker_genes


def explain_model(rf, X_test, marker_genes):
    """Generate SHAP summary plot for Random Forest model."""
    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X_test)
    shap.summary_plot(
        shap_values,
        X_test,
        feature_names=marker_genes,
        plot_type="bar",
        max_display=20,
        show=False,
    )
    plt.savefig(os.path.join(FIGURES_DIR, "shap_summary_plot.png"))
    plt.close()
