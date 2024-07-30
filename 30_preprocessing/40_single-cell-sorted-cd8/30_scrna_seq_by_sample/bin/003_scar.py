#!/usr/bin/env python3

"""
Usage:
 003_scar.py --filtered_adata=<filtered_adata> --id=<id>[options]

Mandatory arguments:
    --filtered_adata=<filtered_adata>        adata filtered
    --id=<id>
        
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

filtered_adata = args["--filtered_adata"]
sample_id = args["--id"]

input_path = "/data/projects/2021/MicrobialMetabolites/single-cell-sorted-cd8/scrna_seq_data/samples" #in the config
sample_id_cleaned = sample_id.rstrip("_adata_filtered")
raw_adata_path = f"{input_path}/adata_{sample_id_cleaned}.h5ad"

print(raw_adata_path)

cpus = 8

# Raw adata 
raw_adata = scvi.data.read_h5ad(raw_adata_path)

#Filter adata
adata = scvi.data.read_h5ad(filtered_adata)

raw_adata.obs["value"] = 0
adata.obs["value"] = 0

scvi.external.SCAR.setup_anndata(adata)

scvi.external.SCAR.get_ambient_profile(
    adata=adata, raw_adata=raw_adata, prob=0.5
)

model = scvi.external.SCAR(adata)
model.train(early_stopping=True, use_gpu=False)

adata.obsm["X_scAR"] = model.get_latent_representation()
adata.layers["denoised"] = model.get_denoised_counts()


# Save raw counts to layer counts and use denoised counts for subsequent analysis
adata.layers["counts"] = adata.X.copy()
adata.X = adata.layers["denoised"]

adata.write_h5ad(f"{sample_id}_denoised_adata.h5ad")