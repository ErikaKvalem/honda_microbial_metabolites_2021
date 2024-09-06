library(readr)
library(ggplot2)
library(ggsignif)
library(forcats)
library(tidyr)
library(ggpubr)
library(rstatix)
library("biomaRt")
conflicts_prefer(stats::filter)
conflicts_prefer(dplyr::rename)



library(dplyr)
# #####################
path="/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/40_tables/40_single-cell-sorted-cd8/40_gex_surface_prot/"
log1p_norm_counts = read.table(paste0(path,"log1p_norm_counts.csv"),sep= ",", header=TRUE, row.names=1)
#library(data.table) 
#log1p_norm_counts <- fread(paste0(path,"log1p_norm_counts.csv"))
samplesheet = read_csv(paste0(path,"samplesheet.csv"))
#samplesheet$sample_id <- gsub("-", "_", samplesheet$sample_id)

counts <- log1p_norm_counts

colnames(counts) <- sub("^X", "", colnames(counts))

counts$gene_id <- rownames(counts)

# Optionally, move "gene_id" to the first column
counts <- counts[, c("gene_id", setdiff(names(counts), "gene_id"))]

resDir = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/40_tables/40_single-cell-sorted-cd8/40_gex_surface_prot"

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ensembl_ids <- c(rownames(counts)) 
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_ids,
                      mart = ensembl)
colnames(gene_symbols)[1] <- "gene_id"
colnames(gene_symbols)[2] <- "gene_name"

# Merge the gene symbols with the original dataframe 'counts'
counts <- merge(counts, gene_symbols, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
                
###################################
# Convert log1p_norm_counts to dataframe
df_nc <- as.data.frame(counts)
library(dplyr)

# Assuming your data frame is called df
df_nc <- df_nc %>%
  select(-gene_name.y) %>%  # Remove the gene_name.y column
  rename(gene_name = gene_name.x)  # Rename gene_name.x to gene_name

# List of gene names
gene_name_value = "Cxcr3"
gene_id_value <- df_nc %>% 
  filter(gene_name.x == gene_name_value) %>% 
  pull(gene_id)



# Initialize an empty dataframe to store the result
result <- data.frame()
# Initialize an empty dataframe to store gene names and p-values
p_values_df <- data.frame(gene_name = character(), p_val = numeric())

df_nc_filtered <- df_nc[df_nc$gene_id == gene_id_value, ]

# Transpose the filtered dataframe
df_t <- t(df_nc_filtered)

# Convert transposed dataframe to data frame
df_t <- as.data.frame(df_t)
###
colnames(df_t)[colnames(df_t) == colnames(df_t)] <- "log1p_norm"
df_t$id <- rownames(df_t)
rownames(df_t) <- NULL
df_t <- df_t[-1, ]
rownames(df_t) <- df_t$id

# Optionally, remove the 'id' column from the dataframe if you no longer need it as a separate column
df_t$id <- NULL
df_t <- subset(df_t, rownames(df_t) != "gene_name")

####
df_t$log1p_norm <- as.numeric(df_t$log1p_norm)

# Add 'gene_name' column
df_t$gene_id <- as.character(gene_id_value)
df_t$gene_name <- as.character(gene_name_value)
df_t$sample_ID<- rownames(df_t)

samplesheet$sample_ID <- gsub("-", "_", samplesheet$index)



df_t <- merge(df_t, samplesheet[, c("sample_ID", "condition", "cell_type")], by = "sample_ID")
# Filter for male
baseline_df <- subset(df_t, condition == "10mix")

# Filter for female
perturbation_df <- subset(df_t, condition == "11mix")


desired_order <- c("10mix", "11mix","GF", "GF-plus")
facet_labels <- c("COLON_Naive","COLON_Intermediate","COLON_Infg","COLON_Exhausted","MPEC_Effector" ,
                  "MPEC_Intermediate",  "MPEC_Progenitor",  "SLEC_Progenitor" , "SLEC_Plastic"
                  , "SLEC_Intermediate", "SLEC_Effector","SLEC_Inf" ,"SLEC_Terminal"  )
# Append the result
result <- rbind(result, df_t)
colnames(result)[2] <- "log1p_norm"
  

result <- result[!is.na(result$condition), ]
result$"condition" <- factor(result$"condition", levels = desired_order)
p  <-  ggviolin(result, x = "condition", y = "log1p_norm", color="condition", trim=FALSE, title=paste0(gene_name_value, " gene expression")) +  
  #facet_wrap(~ cell_type) 
  facet_wrap(~ cell_type, labeller = as_labeller(setNames(facet_labels, unique(result$cell_type)))) 

max_y <- max(result$log1p_norm, na.rm = TRUE)

p <- p + geom_signif(
  comparisons = list(c("10mix", "11mix"), c("11mix", "GF"), c("10mix","GF"),  c("GF-plus","GF")),
  map_signif_level = TRUE,
  test = "wilcox.test",
  y_position = c(5, 6.0,6.5,7.5),  # Different heights for each comparison
  tip_length = 0, 
  vjust = 0.2
) + ylim(0, 8)
  #geom_jitter(shape=16, position=position_jitter(0.2)) + ylim(0, 8) 
p
path_figures = "/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/20_scripts/40_single-cell-sorted-cd8/40_gex_surface_prot/figures/"

ggsave(paste0(path_figures,"gene_expression.png"), p, width =16, height = 14)
