# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

# %%
import pandas as pd
import dorothea
import progeny
import scanpy as sc
import scanpy_helpers as sh
from pathlib import Path
import altair as alt
import numpy as np

sc.set_figure_params(figsize=(4, 4))

# %% [markdown]
# # Build AnnData object

# %%
expr = pd.read_csv(
    "../data/new_data/01_rnaseq_pipeline/star_salmon/salmon.merged.gene_tpm.tsv",
    sep="\t",
)

# %%
adata = sc.AnnData(X=expr.drop(columns=["gene_id"]).set_index("gene_name").T)

# %%
adata.obs["organoid"] = [x[0] for x in adata.obs.index.str.split("_")]
adata.obs["condition"] = [x[1] for x in adata.obs.index.str.split("_")]

# %%
# # flip (likely) switched samples
# adata.obs["condition"] = [
#     {"10mix": "11mix", "11mix": "10mix"}.get(c, c) if o == "mHCO3" else c
#     for o, c in zip(adata.obs["organoid"], adata.obs["condition"])
# ]

# %%
adata.obs

# %% [markdown]
# # PCA

# %%
sc.pp.log1p(adata)

# %% [markdown]
# ### Original

# %%
sc.tl.pca(adata)
sc.pl.pca(adata, color=["organoid", "condition"], size=120, wspace=0.7)

# %% [markdown]
# ### regressed out

# %%
adata_regressed_out = sc.pp.regress_out(adata, "organoid", copy=True)

# %%
sc.tl.pca(adata_regressed_out)

# %%
sc.pl.pca(adata_regressed_out, color=["organoid", "condition"], size=120, wspace=0.7)

# %% [markdown]
# # GSEA results

# %%
de_res_path = Path(
    "/data/projects/2021/MicrobialMetabolites/new_data/deseq2icbi_switched_sample/deseq2_11mix_vs_10mix"
)

# %%
ora_res = []
for pathway_set in ["GO_BP", "GO_MF", "KEGG", "Reactome"]:
    ora_res.append(
        pd.read_csv(
            de_res_path / f"11mix_10mix_ORA_{pathway_set}.tsv", sep="\t"
        ).assign(db=pathway_set)
    )

# %%
ora_res_filtered = (
    pd.concat(ora_res)
    .loc[lambda x: x["p.adjust"] < 0.1]
    .assign(log_fdr=lambda x: -np.log10(x["p.adjust"]))
)

# %%
alt.Chart(ora_res_filtered).mark_circle(size=120, opacity=1).encode(
    x=alt.X(
        "log_fdr", scale=alt.Scale(domain=[1, 2]), axis=alt.Axis(title="log10(FDR)")
    ),
    y=alt.Y("Description", sort="-x", axis=alt.Axis(labelLimit=500, title=None)),
    color="db",
)

# %% [markdown]
# # Progeny

# %%
adata.var_names_make_unique()

# %%
mod = progeny.load_model(organism="Mouse")

# %%
progeny.run(adata, mod)

# %%
ad_progeny = progeny.extract(adata)
sc.pp.regress_out(ad_progeny, "organoid")
ad_progeny = ad_progeny[ad_progeny.obs["condition"].str.contains("mix"), :].copy()

# %%
ad_progeny.obs

# %%
sh.compare_groups.lm.test_lm(
    ad_progeny,
    "~ C(condition, Treatment('10mix')) + organoid",
    groupby="condition",
    contrasts="Treatment('10mix')",
).pipe(sh.util.fdr_correction).sort_values("pvalue")

# %%
heatmap_df = (
    pd.DataFrame(ad_progeny.X, columns=ad_progeny.var_names, index=ad_progeny.obs_names)
    .join(ad_progeny.obs)
    .melt(id_vars=["condition", "organoid"], var_name="pathway", value_name="score")
)

# %%
heatmap_df

# %%
alt.Chart(heatmap_df).mark_rect().encode(
    x="pathway",
    y="organoid",
    color=alt.Color("score", scale=alt.Scale(scheme="redblue", domain=[-1.5, 1.5])),
).facet(row="condition", spacing=5)

# %%
sh.pairwise.plot_paired(
    ad_progeny[ad_progeny.obs["condition"].str.contains("mix"), :],
    var_names=["EGFR", "VEGF", "JAK-STAT", "TGFb"],
    groupby="condition",
    paired_by="organoid",
    pvalue_template=lambda x: "p={:.2f}, paired t-test".format(x),
)

# %% [markdown]
# # Dorothea

# %%
reg = dorothea.load_regulons(organism="Mouse", levels=["A", "B"])

# %%
dorothea.run(adata, reg, scale=True)

# %%
ad_doro = dorothea.extract(adata)
sc.pp.regress_out(ad_doro, "organoid")
ad_doro = ad_doro[ad_doro.obs["condition"].str.contains("mix"), :].copy()

# %%
ad_doro

# %%
doro_lm_res = (
    sh.compare_groups.lm.test_lm(
        ad_doro,
        "~ C(condition, Treatment('10mix')) + organoid",
        groupby="condition",
        contrasts="Treatment('10mix')",
    )
    .dropna()
    .pipe(sh.util.fdr_correction)
    .sort_values("pvalue")
)

# %%
doro_lm_res.head()

# %%
heatmap_df = (
    pd.DataFrame(ad_doro.X, columns=ad_doro.var_names, index=ad_doro.obs_names)
    .join(ad_doro.obs)
    .melt(id_vars=["condition", "organoid"], var_name="TF", value_name="score")
).loc[lambda x: x["TF"].isin(doro_lm_res[:15]["variable"])]

# %%
alt.Chart(heatmap_df).mark_rect().encode(
    x="TF",
    y="organoid",
    color=alt.Color("score", scale=alt.Scale(scheme="redblue")),
).facet(row="condition", spacing=5)

# %%
sh.pairwise.plot_paired(
    ad_doro[ad_doro.obs["condition"].str.contains("mix"), :],
    groupby="condition",
    paired_by="organoid",
    var_names=["Nfkb2", "Rxra", "Foxo3", "Ppara"],
    n_cols=8,
    pvalue_template=lambda x: "p={:.2f}, paired t-test".format(x),
)
