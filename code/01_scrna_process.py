#!/usr/bin/env python3
"""scRNA-seq processing for GBM samples GSM4141788/89/90 (GSE138794).

Input: dense cells x genes count matrices in data/raw/*_dense.csv.gz.
Each file's first column is the cell barcode and the second column ("CLUSTER")
is the original authors' cluster id, which we keep as metadata but do not use
to drive the analysis.

Pipeline: load -> per-cell QC -> filter -> normalize/log1p -> HVG ->
scale/PCA -> neighbors -> Leiden clustering -> UMAP -> marker genes ->
coarse cell-type annotation. Outputs an .h5ad and figures/tables under
results/scrna/.

Run:
    python code/01_scrna_process.py
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "results" / "scrna"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1
sc.settings.figdir = FIG
sc.settings.set_figure_params(dpi=120, frameon=False)

SAMPLES = {
    "GSM4141788_92017": RAW / "GSM4141788_92017_dense.csv.gz",
    "GSM4141789_92217": RAW / "GSM4141789_92217_dense.csv.gz",
    "GSM4141790_92717": RAW / "GSM4141790_92717_dense.csv.gz",
}


def load_sample(name: str, path: pathlib.Path) -> sc.AnnData:
    """Read one dense cells x genes csv.gz into an AnnData."""
    df = pd.read_csv(path, index_col=0)
    orig_cluster = df.pop("CLUSTER").astype(str) if "CLUSTER" in df.columns else None
    # Barcodes are integers in the file; make them unique per sample as strings.
    df.index = [f"{name}_{bc}" for bc in df.index]
    # Counts are sparse (mostly zeros); store as CSR to keep memory/disk small.
    adata = sc.AnnData(
        X=sparse.csr_matrix(df.values.astype(np.float32)),
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns),
    )
    adata.obs["sample"] = name
    if orig_cluster is not None:
        adata.obs["orig_cluster"] = orig_cluster.values
    return adata


def main() -> None:
    adatas = [load_sample(n, p) for n, p in SAMPLES.items()]
    # outer join keeps the union of genes; fill_value=0 is required so genes
    # absent from a sample stay zero counts (the default NaN fill breaks QC).
    adata = sc.concat(adatas, join="outer", label="batch", index_unique=None, fill_value=0)
    adata.obs_names_make_unique()
    print(f"Combined: {adata.n_obs} cells x {adata.n_vars} genes")
    print(adata.obs["sample"].value_counts())

    # --- QC metrics ---
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        groupby="sample",
        jitter=0.4,
        rotation=45,
        multi_panel=True,
        show=False,
        save="_qc_violin.png",
    )

    n_before = adata.n_obs
    # Filter cells: require a reasonable gene count and bounded mito fraction.
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs["n_genes_by_counts"] < 7000].copy()
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    # Filter genes detected in too few cells.
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"QC: {n_before} -> {adata.n_obs} cells, {adata.n_vars} genes")

    # Keep raw counts for downstream DE / reference.
    adata.layers["counts"] = adata.X.copy()

    # --- Normalize ---
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # --- HVG (per-sample batch-aware) ---
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, flavor="seurat", batch_key="sample"
    )
    sc.pl.highly_variable_genes(adata, show=False, save="_hvg.png")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    # --- Scale / PCA ---
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")
    sc.pl.pca_variance_ratio(adata_hvg, n_pcs=50, show=False, save="_pca_variance.png")

    # --- Neighbors / clustering / UMAP ---
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata_hvg, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata_hvg)

    # Carry results back onto the full (all-genes) object.
    adata.obs["leiden"] = adata_hvg.obs["leiden"]
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.uns["leiden"] = adata_hvg.uns.get("leiden", {})

    print("Leiden clusters:")
    print(adata.obs["leiden"].value_counts().sort_index())

    sc.pl.umap(
        adata,
        color=["leiden", "sample", "pct_counts_mt", "n_genes_by_counts"],
        wspace=0.4,
        ncols=2,
        show=False,
        save="_overview.png",
    )

    # --- Marker genes per cluster ---
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=False, save="_markers.png")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(OUT / "cluster_markers.csv", index=False)

    # --- Coarse cell-type annotation by canonical GBM TME markers ---
    marker_sets = {
        "Malignant/Astrocyte": ["EGFR", "GFAP", "SOX2", "OLIG1", "OLIG2", "VIM", "CHI3L1"],
        "Myeloid/Microglia": ["PTPRC", "CD68", "AIF1", "CSF1R", "C1QA", "C1QB", "ITGAM"],
        "T cell": ["CD3D", "CD3E", "CD2", "CD8A"],
        "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MAG"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5"],
        "Pericyte/Mural": ["PDGFRB", "RGS5", "ACTA2"],
    }
    present = {k: [g for g in v if g in adata.raw.var_names] for k, v in marker_sets.items()}
    for ct, genes in present.items():
        if genes:
            sc.tl.score_genes(adata, genes, score_name=f"score_{ct}", use_raw=True)

    score_cols = [f"score_{ct}" for ct in present if present[ct]]
    scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    scores.columns = [c.replace("score_", "") for c in scores.columns]
    cluster_label = scores.idxmax(axis=1)
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_label).astype("category")
    scores.assign(assigned=cluster_label).to_csv(OUT / "cluster_celltype_scores.csv")

    print("Cluster -> cell type:")
    print(cluster_label)

    sc.pl.umap(adata, color=["leiden", "cell_type"], wspace=0.4, show=False, save="_celltype.png")

    flat_markers = [g for genes in present.values() for g in genes]
    sc.pl.dotplot(adata, flat_markers, groupby="cell_type", show=False, save="_celltype_markers.png", use_raw=True)

    # --- Save ---
    h5ad = OUT / "gbm_scrna_processed.h5ad"
    adata.write(h5ad)
    print(f"Wrote {h5ad}")

    # Summary table.
    summary = (
        adata.obs.groupby(["cell_type", "sample"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    summary.to_csv(OUT / "celltype_by_sample.csv")
    print(summary)


if __name__ == "__main__":
    main()
