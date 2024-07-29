#!/usr/bin/env python3

"""
Usage:
  004_solo.py --denoised_filtered_adata=<denoised_filtered_adata>   [options]

Mandatory arguments:
    --denoised_filtered_adata=<denoised_filtered_adata>                raw
 
        
Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import matplotlib.pyplot as plt

import scvi
import scanpy as sc
import pandas as pd


args = docopt(__doc__)
denoised_filtered_adata_path = args["--denoised_filtered_adata"]


adata_denoised_filtered = scvi.data.read_h5ad(denoised_filtered_adata_path)

sc.pp.highly_variable_genes(
    adata_denoised_filtered,
    n_top_genes=4000,
    subset=False,
    flavor="seurat_v3",
    batch_key="sample_id"
)

adata_scvi = adata_denoised_filtered[:, adata_denoised_filtered.var["highly_variable"]].copy()


scvi.model.SCVI.setup_anndata(
    adata_scvi,
    layer="denoised",
   batch_key="sample_id"
)

model = scvi.model.SCVI(adata_scvi) #set up scvi model for later on predict doublets 
model.train(early_stopping=True, use_gpu=True)

adata_raw = adata_denoised_filtered.copy()
sc.pp.normalize_total(adata_raw)
sc.pp.log1p(adata_raw)
sc.tl.pca(adata_raw, use_highly_variable=True)
adata_denoised_filtered.raw =adata_raw

adata_denoised_filtered.obsm["X_pca"] = adata_raw.obsm["X_pca"]

sc.pp.neighbors(adata_denoised_filtered, key_added="neighbors_uncorrected")
sc.tl.umap(adata_denoised_filtered, neighbors_key="neighbors_uncorrected")
adata_denoised_filtered.obsm["X_umap_uncorrected"] = adata_denoised_filtered.obsm["X_umap"]
del adata_denoised_filtered.obsm["X_umap"]

adata_denoised_filtered.obs


sc.pl.embedding(adata_denoised_filtered, "umap_uncorrected", color=["sample_id", "n_genes", "pct_counts_mt"], cmap="inferno")


# ## With batch correction
adata_denoised_filtered.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(adata_denoised_filtered, use_rep="X_scVI")
sc.tl.umap(adata_denoised_filtered)

sc.pl.umap(adata_denoised_filtered, color=["sample_id", "n_genes", "pct_counts_mt"], cmap="inferno")

# ## Doublet detection (SOLO)
# train solo model for doublet detection
solo_models = [
    scvi.external.SOLO.from_scvi_model(model, restrict_to_batch=batch)
    for batch in adata_denoised_filtered.obs["sample_id"].unique()
]

def run_solo(solo_model):
    solo_model.train()
    return solo_model.predict(soft=False)

solo_res = []
for solo_model in solo_models:
    solo_res.append(run_solo(solo_model))

adata_denoised_filtered.obs["is_doublet"] = pd.concat(solo_res) #predicted label is_doublet

solo_series = pd.concat(solo_res)

solo_series.index = solo_series.index.str.replace("-0$", "", regex=True) #remove dash 0 that scvi adds to each barcode 

adata_denoised_filtered.obs["is_doublet"] = solo_series 

#Is doublet umap 
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sc.pl.umap(adata_denoised_filtered, color="is_doublet", ax=ax)
fig.savefig("is_doublet.png", bbox_inches="tight")

#subset adata to keep only barcodes that are label singlet 
adata_nodoublet = adata_denoised_filtered[adata_denoised_filtered.obs["is_doublet"] == "singlet", :].copy()

#recalculate the scvi model once the doublet have been filtered
adata_nodoublet.obsm["X_scVI"] = model.get_latent_representation(adata_scvi[adata_nodoublet.obs_names, :])

sc.pp.neighbors(adata_nodoublet, use_rep="X_scVI")

sc.tl.umap(adata_nodoublet)

sc.tl.leiden(adata_nodoublet, key_added="leiden", resolution=0.3)

# Save adata nodoublet
adata_nodoublet.write_h5ad("adata_nodoublet.h5ad")