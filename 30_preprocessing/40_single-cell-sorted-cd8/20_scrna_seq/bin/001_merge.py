#!/usr/bin/env python3

"""
Usage:
  001_merge.py --count_dir=<count_dir> --experiment=<experiment> [options]

Mandatory arguments:
    --count_dir=<count_dir>                count directory
    --experiment=<experiment>
         

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
from functools import reduce
from pathlib import Path
import scipy

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from tqdm.contrib.concurrent import process_map

args = docopt(__doc__)
count_dir = args["--count_dir"]
experiment = args["--experiment"]

if experiment == "2019-10-29_sorted_cd8":
    sample_names = ["10mix1","10mix2","11mix1","11mix2","GF1","GF2"]
    
elif experiment =="2021-02-01_sorted_cd8_til":
    sample_names = ["10mix-ICI1","10mix-ICI2","11mix-ICI1","11mix-ICI2","GF-ICI1","GF-ICI2","GF-ICI1-plus","GF-ICI2-plus"]


# Loop through each sample name and create empty DataFrames
for sample_name in sample_names:

        
    path_barcodes = f"{count_dir}/barcodes.tsv.gz"
    path_features = f"{count_dir}/features.tsv.gz"
    path_matrix = f"{count_dir}/matrix.mtx.gz"
    
    # Create an empty DataFrame for df_barcodes
    df_barcodes = pd.read_csv(path_barcodes, delimiter="\t", header=None)
    df_barcodes = df_barcodes.rename(columns= {0:"barcodes"})
    df_barcodes.index = df_barcodes["barcodes"].str.split("-").str[0]
    df_barcodes = df_barcodes.drop(["barcodes"],axis=1)
    df_barcodes.index.name = None
    # Create an empty DataFrame for df_features
    df_features = features = pd.read_csv(path_features, delimiter="\t", header=None, compression = "gzip", index_col=1)
    df_features.index.name = None
    df_features = df_features.rename(columns={0:"ensembl_id",2:"feature_types"})
    
    mat = scipy.io.mmread(path_matrix)
    adata = sc.AnnData(X=mat.T)
    adata.X = scipy.sparse.csr_matrix(adata.X)
    
    adata.var = pd.DataFrame(df_features)
    adata.var["gene_name"]=adata.var_names
    adata.var_names = adata.var["ensembl_id"]
    adata.var.index.name = None
    assert  adata.var_names.is_unique
    adata.obs = pd.DataFrame(df_barcodes, index=df_barcodes.index.values)
    adata.obs["sample_id"]=sample_name
    assert  adata.obs_names.is_unique
    
    # Save the AnnData object to a h5ad file
    adata.write(f'adata_{sample_name}.h5ad')


adata_list = []

# Loop through each sample name and create empty DataFrames
for sample_name in sample_names:
    file_path = f'adata_{sample_name}.h5ad'
    adata = ad.read_h5ad(file_path)
    adata_list.append(adata)
    
# Concatenate the AnnData objects along the rows (observations)
merged_adata = adata_list[0].concatenate(adata_list[1:], join="outer")
merged_adata.X = csr_matrix(merged_adata.X)
merged_adata.write(f'{experiment}_merged_data.h5ad')
