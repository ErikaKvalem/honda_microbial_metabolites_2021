########### This code integrates scFEA outputs and Seurat object of the neutrophil analysis ##########
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)

ad_path <- "adata_solo_annotated_all.h5ad"
sceasy::convertFormat(ad_path, from="anndata", to="seurat", outFile="adata_solo_annotated_all.rds")

##### load Seurat object 
obj <- readRDS(file = "adata_solo_annotated_all.rds")
df_flux <- read.csv(paste0("output/mouse_flux.csv"),header = TRUE, row.names = 1)
df_balance <- read.csv(paste0("output/mouse_balance.csv"),header = TRUE, row.names = 1)



##### cluster based on fluxome counts instead of RNA counts 
# add flux and balance as a new assay
predFlux <- data.matrix(df_flux)
predFlux0 <- t(predFlux)

balance <- data.matrix(df_balance)
balance0 <- t(balance)

obj[["FLUX"]] <- CreateAssayObject(counts = predFlux0)
obj[["BALANCE"]] <- CreateAssayObject(counts = balance0)

# cluster based on FLUX by defining FLUX as default assay instead of RNA
DefaultAssay(obj) <- 'FLUX'
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = F)
obj <- ScaleData(obj, features = rownames(obj), assay = 'FLUX', verbose = F)
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 10, reduction.name = 'pca.flux', verbose = F)
ElbowPlot(obj)
obj <- FindNeighbors(obj, dims = 1:10, verbose = F)
obj <- FindClusters(obj, resolution = 0.2, verbose = F, graph.name = "RNA_nn")
obj <- RunUMAP(obj, dims = 1:10, assay = 'FLUX', reduction.name = "umap.flux", verbose = F)
#these are the flux clusters 
DimPlot(obj, reduction = "umap.flux", group.by = "seurat_clusters", label = TRUE,label.size = 10) 

#overlay flux clusters with umap generated using "RNA" gene counts as DefaultAssay 
DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 10) 

##### find modules and balance compounds that define the flux clusters 
#therefore I'm doing DEG analysis as I usually do for RNA gene counts but for FLUX and BALANCE to find differentially expressed Modules and compounds 
Idents(obj) <- "seurat_clusters"
DefaultAssay(obj) <- "FLUX"
#use RNA assay and data which is the log noramlized counts 
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))
write.csv(markers, "output/DEG_modules_per_cluster.csv")

DefaultAssay(obj) <- "BALANCE"
#use RNA assay and data which is the log noramlized counts 
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))
write.csv(markers, "output/DEG_balance_per_cluster.csv")

## plot some modules and compounds that are differentially expressed between clusters in a feature plot 

#bm cluster 2 specific 
DefaultAssay(obj) <- "FLUX"
FeaturePlot(obj, features = "M-30", reduction = "umap")

#I then checked which genes define module M-61 and doublechecked if the signal of some of the enzyme genes are also clustered in the same umap cluster
DefaultAssay(obj) <- "RNA"
FeaturePlot(obj, features = "UDP.N.acetylglucosamine", reduction = "umap")


###################
library(patchwork)
p1 <- DimPlot(obj, reduction = "umap", group.by = "sample_id") + ggtitle('umap of Gene')
p2 <- DimPlot(obj, reduction = "umap.flux", group.by = "sample_id") + ggtitle('umap of Flux')
p2 + p1

######################
# UMAP plot colored by sample_id
p1 <- DimPlot(obj, reduction = "umap", group.by = "sample_id") + 
  ggtitle("UMAP of Gene by Sample ID")

# UMAP plot colored by cell_type
p2 <- DimPlot(obj, reduction = "umap", group.by = "cell_type") + 
  ggtitle("UMAP of Gene by Cell Type")
p2 + p1


