#!/usr/bin/env python3

"""
Usage:
  002_scar.py --raw_adata=<raw_adata>  --filtered_adata=<filtered_adata> [options]

Mandatory arguments:
    --raw_adata=<raw_adata>                raw
    --filtered_adata=<filtered_adata>        adata filtered
        
Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""


import os
from docopt import docopt
from functools import reduce
import anndata as ad
import pandas as pd
import scvi


args = docopt(__doc__)
raw_adata_path = args["--raw_adata"]
filtered_adata_path = args["--filtered_adata"]
#experiment = args["--experiment"]

cpus = 8

# Raw adata 
raw_adata = scvi.data.read_h5ad(raw_adata_path)

#Filter adata
adata = scvi.data.read_h5ad(filtered_adata_path)


adatas = dict()
for s in adata.obs["sample_id"].unique():

    raw_adata_s = raw_adata[raw_adata.obs["sample_id"] == s, :].copy()
    raw_adata_s.obs["value"] = 0

    adata_s = adata[adata.obs["sample_id"] == s, :].copy()
    adata_s.obs["value"] = 0

    scvi.external.SCAR.setup_anndata(adata_s)

    scvi.external.SCAR.get_ambient_profile(
        adata=adata_s, raw_adata=raw_adata_s, prob=0.8
    )
 
    model = scvi.external.SCAR(adata_s)
    model.train(early_stopping=True, use_gpu=True)

    adata_s.obsm["X_scAR"] = model.get_latent_representation()
    adata_s.layers["denoised"] = model.get_denoised_counts()

    adatas[s] = adata_s

# Integrate back
adata = ad.concat(adatas, join="outer")

# Save raw counts to layer counts and use denoised counts for subsequent analysis
adata.layers["counts"] = adata.X.copy()
adata.X = adata.layers["denoised"]

adata.write_h5ad(f"denoised_adata.h5ad")
