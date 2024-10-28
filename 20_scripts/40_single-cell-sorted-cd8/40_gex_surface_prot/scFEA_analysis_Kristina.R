########### This code integrates scFEA outputs and Seurat object of the neutrophil analysis ##########

##### load Seurat object 
obj <- readRDS(file = "/data/khandl/eos_tumor/seurat_objects/tumor_bm_blood_neutrophils.rds")

##### read in scFEA output files (I run it separate for each cluster)
# flux 
df0 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",0,"_flux.csv"),header = TRUE, row.names = 1)
df1 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",1,"_flux.csv"),header = TRUE, row.names = 1)
df2 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",2,"_flux.csv"),header = TRUE, row.names = 1)
df3 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",3,"_flux.csv"),header = TRUE, row.names = 1)
df4 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",4,"_flux.csv"),header = TRUE, row.names = 1)
df5 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",5,"_flux.csv"),header = TRUE, row.names = 1)
df6 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",6,"_flux.csv"),header = TRUE, row.names = 1)
df7 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",7,"_flux.csv"),header = TRUE, row.names = 1)
df8 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",8,"_flux.csv"),header = TRUE, row.names = 1)
df9 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",9,"_flux.csv"),header = TRUE, row.names = 1)
df10 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",10,"_flux.csv"),header = TRUE, row.names = 1)

#combine to one dataframe 
df_flux <- rbind(df0,df1)
df_flux <- rbind(df_flux,df2)
df_flux <- rbind(df_flux,df3)
df_flux <- rbind(df_flux,df4)
df_flux <- rbind(df_flux,df5)
df_flux <- rbind(df_flux,df6)
df_flux <- rbind(df_flux,df7)
df_flux <- rbind(df_flux,df8)
df_flux <- rbind(df_flux,df9)
df_flux <- rbind(df_flux,df10)

# balance 
df0 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",0,"_balance.csv"),header = TRUE, row.names = 1)
df1 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",1,"_balance.csv"),header = TRUE, row.names = 1)
df2 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",2,"_balance.csv"),header = TRUE, row.names = 1)
df3 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",3,"_balance.csv"),header = TRUE, row.names = 1)
df4 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",4,"_balance.csv"),header = TRUE, row.names = 1)
df5 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",5,"_balance.csv"),header = TRUE, row.names = 1)
df6 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",6,"_balance.csv"),header = TRUE, row.names = 1)
df7 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",7,"_balance.csv"),header = TRUE, row.names = 1)
df8 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",8,"_balance.csv"),header = TRUE, row.names = 1)
df9 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",9,"_balance.csv"),header = TRUE, row.names = 1)
df10 <- read.csv(paste0("/data/khandl/scFEA/neut_cluster",10,"_balance.csv"),header = TRUE, row.names = 1)

df_balance <- rbind(df0,df1)
df_balance <- rbind(df_balance,df2)
df_balance <- rbind(df_balance,df3)
df_balance <- rbind(df_balance,df4)
df_balance <- rbind(df_balance,df5)
df_balance <- rbind(df_balance,df6)
df_balance <- rbind(df_balance,df7)
df_balance <- rbind(df_balance,df8)
df_balance <- rbind(df_balance,df9)
df_balance <- rbind(df_balance,df10)


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
write.csv(markers, "/scratch/khandl/eos_tumor/neutrophils/scFEA/DEG_modules_per_cluster.csv")

DefaultAssay(obj) <- "BALANCE"
#use RNA assay and data which is the log noramlized counts 
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))
write.csv(markers, "/scratch/khandl/eos_tumor/neutrophils/scFEA/DEG_balance_per_cluster.csv")

## plot some modules and compounds that are differentially expressed between clusters in a feature plot 

#bm cluster 2 specific 
DefaultAssay(obj) <- "FLUX"
FeaturePlot(obj, features = "M-61", reduction = "umap")

#I then checked which genes define module M-61 and doublechecked if the signal of some of the enzyme genes are also clustered in the same umap cluster
DefaultAssay(obj) <- "RNA"
FeaturePlot(obj, features = "Aldh2", reduction = "umap")

# blood cluster 1 specific 
DefaultAssay(obj) <- "FLUX"
FeaturePlot(obj, features = "M-139", reduction = "umap")

DefaultAssay(obj) <- "RNA"
FeaturePlot(obj, features = "Dck", reduction = "umap")

# colon/tumor cluster 0 specific 
DefaultAssay(obj) <- "FLUX"
FeaturePlot(obj, features = "M-155", reduction = "umap")

DefaultAssay(obj) <- "RNA"
FeaturePlot(obj, features = "Entpd1", reduction = "umap")

# BM cluster 3 specific  
DefaultAssay(obj) <- "FLUX"
FeaturePlot(obj, features = "M-109", reduction = "umap")

DefaultAssay(obj) <- "RNA"
FeaturePlot(obj, features = "Ugdh", reduction = "umap")

