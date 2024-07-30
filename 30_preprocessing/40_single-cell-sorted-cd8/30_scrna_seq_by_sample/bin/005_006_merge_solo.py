#!/usr/bin/env python3

"""
Usage:
 005_006_merge_solo.py --path_adata_denoised_after_qc=<path_adata_denoised_after_qc> [options]

Mandatory arguments:
  --path_adata_denoised_after_qc=<path_adata_denoised_after_qc>                 path

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import os
import scanpy as sc
import scvi
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad
#import scanpy_helpers as sh
#import decoupler as dc
#from tqdm import tqdm
#import re 


########################## MERGE 
args = docopt(__doc__)
path_adata_denoised_after_qc =args["--path_adata_denoised_after_qc"]

#list of all files in the directory
file_list = os.listdir(path_adata_denoised_after_qc)

adata_list = []
for adata_file in file_list:
    adata = sc.read_h5ad(f"{path_adata_denoised_after_qc}/{adata_file}")
    adata_list.append(adata)

    # Concatenate the AnnData objects along the rows (observations)
merged_adata = ad.concat(adata_list, index_unique="-", join="outer", fill_value=0)

merged_adata.write_h5ad("merged_adata.h5ad")

############################# SOLO

sc.pp.highly_variable_genes(
    merged_adata,
    layer="counts",
    n_top_genes=4000,
    subset=False,
    flavor="seurat_v3",
    batch_key="sample_id"
)


adata_scvi = merged_adata[:, merged_adata.var["highly_variable"]].copy()




#WARNING  No GPU/TPU found, falling back to CPU. (Set TF_CPP_MIN_LOG_LEVEL=0 and rerun for more info.)
scvi.model.SCVI.setup_anndata(
    adata_scvi,
    layer="denoised",
   batch_key="sample_id"
)



model = scvi.model.SCVI(adata_scvi) #set up scvi model for later on predict doublets 



model.train(early_stopping=True) 

adata_raw = merged_adata.copy()

sc.pp.normalize_total(adata_raw)
sc.pp.log1p(adata_raw)
sc.tl.pca(adata_raw, use_highly_variable=True)
merged_adata.raw =adata_raw


merged_adata.obsm["X_pca"] = adata_raw.obsm["X_pca"]

sc.pp.neighbors(merged_adata, key_added="neighbors_uncorrected")
sc.tl.umap(merged_adata, neighbors_key="neighbors_uncorrected")
merged_adata.obsm["X_umap_uncorrected"] = merged_adata.obsm["X_umap"]
del merged_adata.obsm["X_umap"]

merged_adata.obs


sc.pl.embedding(merged_adata, "umap_uncorrected", color=["sample_id", "n_genes", "pct_counts_mt"], cmap="inferno")


# ## With batch correction
merged_adata.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(merged_adata, use_rep="X_scVI")
sc.tl.umap(merged_adata)

sc.pl.umap(merged_adata, color=["sample_id", "n_genes", "pct_counts_mt"], cmap="inferno")

# ## Doublet detection (SOLO)
# train solo model for doublet detection
solo_models = [
    scvi.external.SOLO.from_scvi_model(model, restrict_to_batch=batch)
    for batch in merged_adata.obs["sample_id"].unique()
]

def run_solo(solo_model):
    solo_model.train()
    return solo_model.predict(soft=False)

solo_res = []
for solo_model in solo_models:
    solo_res.append(run_solo(solo_model))

merged_adata.obs["is_doublet"] = pd.concat(solo_res) #predicted label is_doublet

solo_series = pd.concat(solo_res)

solo_series.index = solo_series.index.str.replace("-0$", "", regex=True) #remove dash 0 that scvi adds to each barcode 

merged_adata.obs["is_doublet"] = solo_series 

#Is doublet umap 
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sc.pl.umap(merged_adata, color="is_doublet", ax=ax)
fig.savefig("is_doublet.png", bbox_inches="tight")

#subset adata to keep only barcodes that are label singlet 
adata_nodoublet = merged_adata[merged_adata.obs["is_doublet"] == "singlet", :].copy()

#recalculate the scvi model once the doublet have been filtered
adata_nodoublet.obsm["X_scVI"] = model.get_latent_representation(adata_scvi[adata_nodoublet.obs_names, :])

sc.pp.neighbors(adata_nodoublet, use_rep="X_scVI")

sc.tl.umap(adata_nodoublet)

sc.tl.leiden(adata_nodoublet, key_added="leiden", resolution=0.3)

# Save adata nodoublet
adata_nodoublet.write_h5ad("adata_nodoublet.h5ad")