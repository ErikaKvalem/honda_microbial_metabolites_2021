########### This code integrates scFEA outputs and Seurat object of the neutrophil analysis ##########
# Enter commands in R (or R studio, if installed)
#install.packages('Seurat')
library(Seurat)

ad_path <- "adata_solo_re_annotated_all.h5ad"
sceasy::convertFormat(ad_path, from="anndata", to="seurat", outFile="adata_solo_re_annotated_all.rds")

##### load Seurat object 
##### load Seurat object 
obj <- readRDS(file = "adata_solo_re_annotated_all.rds")
df_flux <- read.csv(paste0("output/output_annotated_til_colon/mouse_flux.csv"),header = TRUE, row.names = 1)
df_balance <- read.csv(paste0("output/output_annotated_til_colon/mouse_balance.csv"),header = TRUE, row.names = 1)



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
write.csv(markers, "output/output_re_annotated_til_colon/DEG_modules_per_cluster.csv")


DefaultAssay(obj) <- "BALANCE"
#use RNA assay and data which is the log noramlized counts 
markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC))
write.csv(markers, "output/output_re_annotated_til_colon/DEG_balance_per_cluster.csv")

## plot some modules and compounds that are differentially expressed between clusters in a feature plot 

#bm cluster 2 specific 
DefaultAssay(obj) <- "FLUX"
FeaturePlot(obj, features = "M-30", reduction = "umap")

#I then checked which genes define module M-61 and doublechecked if the signal of some of the enzyme genes are also clustered in the same umap cluster
DefaultAssay(obj) <- "RNA"
FeaturePlot(obj, features = markers$gene, reduction = "umap")


###################
library(patchwork)
p1 <- DimPlot(obj, reduction = "umap", group.by = "cell_type") + ggtitle('umap of Gene')
p2 <- DimPlot(obj, reduction = "umap.flux", group.by = "leiden_res0_5", label=FALSE) + ggtitle('umap of Flux')
p2_labeled <- LabelClusters(p2, id = "leiden_res0_5")

p3 <- DimPlot(obj, reduction = "umap", group.by = "leiden_res0_5", label=FALSE) + ggtitle('umap of Flux')
p3_labeled <- LabelClusters(p3, id = "leiden_res0_5")
p1 + p2_labeled + p3_labeled
###################


library(Seurat)
library(ggpubr)
library(ggplot2)
library(gridExtra)

# Get the expression matrix for all genes
expr_matrix <- FetchData(obj, vars = markers$gene)

# Determine the global min and max values across all genes
global_min <- min(expr_matrix, na.rm = TRUE)
global_max <- max(expr_matrix, na.rm = TRUE)

# Start a PDF device to save the plots
pdf("FeaturePlots.pdf", width = 12, height = 8)

# Create a list to hold the plots
plots <- list()

# Loop over each gene in markers$gene and generate a FeaturePlot with consistent scale
for (gene in markers$gene) {
  p <- FeaturePlot(obj, features = gene, reduction = "umap") + 
    scale_color_gradient(limits = c(global_min, global_max)) + # Apply the global scale
    ggtitle(paste("Feature Plot for", gene))
  
  plots[[gene]] <- p
}

# Arrange plots in a grid with 2 plots per row and multiple pages if needed
grid_plots <- marrangeGrob(grobs = plots, nrow = 2, ncol = 2)

# Save the arranged plots to the PDF
print(grid_plots)

# Close the PDF device
dev.off()

####################
DimPlot(obj, group.by = "cell_type", reduction = "umap", label = FALSE, 
        cols= c("#fac720", "#000000", "#000000", "#000000", "#000000","#000000", "#000000", "#000000", "#000000", "#000000", "#000000", "#000000"))
DimPlot(obj, group.by = "leiden_res0_5", reduction = "umap", label = FALSE, 
        cols= c("#e30800", "#f56505", "#dec400", "#006630", "#0223c7","#5b02c7", "#00b0e6", "#c40080", "#02f00a", "#7d3301", "#000000"))



df <- data.frame(obj@meta.data$origin, obj@meta.data$condition) 
colnames(df) <- c("Cell_Type", "Condition")
df <- df %>% group_by(Condition, Cell_Type) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(Percent = Nb/C*100) 

library("viridis")

df<-df[df$Cell_Type!="Unknown",]
# df$Condition <- gsub("Control", "CTL", df$Condition)
df$Cell_Type <- gsub("Stellates_Mesenchymal", "Stellates", df$Cell_Type)
df$Cell_Type <- gsub("PP_Gamma", "PP", df$Cell_Type)

xtheme <- theme_bw()+ theme(plot.title = element_text(face = "bold" ,hjust = 0.5, size= 10)
                            ,axis.text.y = element_text(face = "bold",angle = 0, size = 10, hjust=1)
                            ,axis.title.y = element_text(face = "bold", size = rel(0.8))
                            ,axis.text.x = element_text(face = "bold",angle = 0, size = 10)
                            ,axis.title.x = element_text(face = "bold", size = rel(0.8))
                            ,axis.ticks.x=element_blank(), strip.text = element_text(size=10)) 
ggplot(df, aes(fill=Condition, y=Percent, x=Cell_Type)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_viridis(discrete = T) + xlab("") + xtheme + 
  theme(legend.position='top', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

##################################


## Median Gene Per Sample
df_mgps <- data.frame(obj@meta.data$nFeaturess_RNA, obj@meta.data$sample_id) 
colnames(df_mgps) <- c("Median_Gene_Number", "Condition")
ggplot(df_mgps, aes(Condition, Median_Gene_Number)) + geom_boxplot(aes(fill = Condition), 
                                                                   width=0.5, outlier.size = 0.2) +
  scale_fill_viridis(discrete = TRUE) + xlab("") + theme(legend.position = "none")

VlnPlot(obj, features = c("nFeaturess_RNA", "nCounts_RNA", "pct_counts_mt"), ncol = 3, group.by = "sample_id",
        pt.size = 0, combine = T) & theme(legend.position = 'none', 
                                          axis.title.x = element_blank())
##################################

library(scSorter)

# Define the gene markers and their corresponding cell types
gene_markers <- c(
  "Ccr7", "Sell", "Tcf7", "Lef1", "Il7r",          # Naive CD8+ T cells
  "Gzmb", "Prf1", "Ifng", "Tnf", "Klrg1", "Cd69",  # Effector CD8+ T cells
  "Pdcd1", "Lag3", "Tigit", "Havcr2", "Eomes", "Tox", "Ctla4",  # Exhausted CD8+ T cells
  "Cd44", "Cd45rb", "Il7r", "Ccr7", "Cx3cr1",      # Memory CD8+ T cells
  "Itgae", "Cd69", "Cxcr6", "Hobit", "Zeb2",       # Tissue-Resident Memory T cells
  "Tcf7", "Lef1", "Cd28", "Cxcr5",                 # Stem-Like CD8+ T cells
  "Foxp3", "Il2ra", "Ctla4", "Entpd1", "Nt5e"      # Regulatory CD8+ T cells
)

cell_types <- c(
  rep("Naive CD8+ T cells", 5),
  rep("Effector CD8+ T cells", 6),
  rep("Exhausted CD8+ T cells", 7),
  rep("Memory CD8+ T cells", 5),
  rep("Tissue-Resident Memory T cells", 5),
  rep("Stem-Like CD8+ T cells", 4),
  rep("Regulatory CD8+ T cells", 5)
)
# Use biomaRt to get the ENSEMBL IDs
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # Mouse dataset

# Query biomaRt for ENSEMBL IDs corresponding to the gene symbols
gene_info <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = anno_df$Marker,
  mart = ensembl
)

# Merge the ENSEMBL IDs with the original data frame
anno_df <- merge(anno_df, gene_info, by.x = "Marker", by.y = "external_gene_name", all.x = TRUE)

# Rename the ensembl_gene_id column to "ensembl"
colnames(anno_df)[colnames(anno_df) == "ensembl_gene_id"] <- "ensembl"



# Write the data frame to a CSV file
write.csv(anno_df, "anno.csv", row.names = FALSE)

# Show the data frame
print(anno_df)


anno <- read.csv("anno.csv")
anno <- na.omit(anno)
expr <- obj@assays$RNA@data
topgenes <- head(obj@assays[["RNA"]]@meta.features[["ensembl_id"]], 3000)
topgene_filter = rowSums(expr[topgenes, ]!=0) > ncol(expr)*.1
topgenes = topgenes[topgene_filter]

## At last, we subset the preprocessed expression data and run scSorter.
picked_genes = unique(c(anno$Marker, topgenes))
expr = expr[rownames(expr) %in% picked_genes, ]

rts <- scSorter(expr, anno)
obj$cell_type <- rts$Pred_Type
