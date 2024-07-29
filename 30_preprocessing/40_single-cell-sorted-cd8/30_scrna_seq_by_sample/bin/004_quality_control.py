#!/usr/bin/env python3

"""
Usage:
  004_quality_control.py --adata_denoised=<adata_denoised> --id=<id>[options]

Mandatory arguments:
    --adata_denoised=<adata_denoised>
    --id=<id>

        
Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import scanpy as sc
import seaborn as sns
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation

args = docopt(__doc__)
adata_denoised = args["--adata_denoised"]
id = args["--id"]
#id =  str(adata_denoised).split("_")[0]


adata_denoised= sc.read_h5ad(adata_denoised)


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


adata_denoised.obs["outlier"] = (
    is_outlier(adata_denoised, "total_counts", 5)
    | is_outlier(adata_denoised, "n_genes_by_counts", 5)
    | is_outlier(adata_denoised, "pct_counts_in_top_20_genes", 5)
)


#adata_denoised.obs.outlier.value_counts()

adata_denoised.obs["mt_outlier"] = is_outlier(adata_denoised, "pct_counts_mt", 3) | (
    adata_denoised.obs["pct_counts_mt"] > 8
)
############################### not sure about this part need to review 
#adata_denoised.obs.mt_outlier.value_counts()
#print(f"Total number of cells: {adata_denoised.n_obs}")
#adata_denoised_merged = adata_denoised[(~adata_denoised.obs.outlier) & (~adata_denoised.obs.mt_outlier)].copy()
#print(f"Number of cells after filtering of low quality cells: {adata_denoised.n_obs}")
###############################


# Save adata_denoised nodoublet

adata_denoised.write_h5ad(f"{id}_after_qc.h5ad")


