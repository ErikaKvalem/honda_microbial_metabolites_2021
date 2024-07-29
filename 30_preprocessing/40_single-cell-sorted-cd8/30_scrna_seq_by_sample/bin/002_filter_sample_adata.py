#!/usr/bin/env python3

"""
Usage:
 002_filter_sample_adata.py --adata=<adata> --sample_id=<sample_id>[options]

Mandatory arguments:
  --adata=<adata>                 adata
  --id=<id>                 

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import pandas as pd
#import scanpy_helpers as sh
import decoupler as dc
from tqdm import tqdm
import re 

args = docopt(__doc__)
adata = args["--adata"]
sample_id = args["--sample_id"]

#name = str(adata)
#pattern = r'adata_(.+?)\.h5ad'
#match = re.search(pattern, name)
#sample_name = match.group(1)

adata = sc.read_h5ad(adata)

# mitochondrial genes
# mitochondrial genes
adata.var["mt"] = adata.var.gene_name.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var.gene_name.str.startswith(("Rps", "Rpl"))
# hemoglobin genes.
adata.var["hb"] = adata.var.gene_name.str.contains("^Hb[^(P)]")


sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=20)

adata.write_h5ad(f"{sample_id}_adata_filtered.h5ad")