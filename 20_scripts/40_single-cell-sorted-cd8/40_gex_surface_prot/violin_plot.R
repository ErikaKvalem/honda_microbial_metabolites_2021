library(readr)
library(ggplot2)
library(ggsignif)
library(forcats)
library(tidyr)
library(ggpubr)
library(rstatix)
library("biomaRt")



library(dplyr)
# #####################
path="/data/scratch/kvalem/projects/2021/honda_microbial_metabolites_2021/40_tables/40_single-cell-sorted-cd8/40_gex_surface_prot/"
log1p_norm_counts = read.table(paste0(path,"log1p_norm_counts_ps.csv"),sep= ",", header=TRUE, row.names=1)
samplesheet = read_csv(paste0(path,"samplesheet_ps.csv"))
#counts = read.table(paste0(path,"counts_ps.csv"),sep= ",", header=TRUE, row.names=1)
counts <- log1p_norm_counts
#colnames(counts) <- sub("^.", "", colnames(counts))
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

# List of gene names

index_list_name <- c("Ifng","Gnai3")

# Initialize an empty dataframe to store the result
result <- data.frame()
# Initialize an empty dataframe to store gene names and p-values
p_values_df <- data.frame(gene_name = character(), p_val = numeric())

# Iterate over each gene name in index_list
for (i in seq_along(index_list_name)) {
  # Filter dataframe based on gene name
  i=1
  df_nc_filtered <- df_nc[df_nc$gene_name == index_list_name[i], ]
  
  # Transpose the filtered dataframe
  df_t <- t(df_nc_filtered)
  
  # Convert transposed dataframe to data frame
  df_t <- as.data.frame(df_t)
  ###
  colnames(df_t)[colnames(df_t) == "1"] <- "log1p_norm"
  df_t$id <- rownames(df_t)
  rownames(df_t) <- NULL
  df_t <- df_t[-1, ]
  rownames(df_t) <- df_t$id
  
  # Optionally, remove the 'id' column from the dataframe if you no longer need it as a separate column
  df_t$id <- NULL
  df_t <- subset(df_t, rownames(df_t) != "gene_name")
  
  ####
  
  # Rename the column to 'log1p_norm'
 # colnames(df_t)[colnames(df_t) == index_list[i]] <- 'log1p_norm'
  df_t$log1p_norm <- as.numeric(df_t$log1p_norm)
  
  df_t$log1p_norm <- scale(df_t$log1p_norm, center = TRUE, scale = TRUE)
  
  # Apply the log(x + 1) transformation
  #df_t$log1p_norm <- log1p(df_t$log1p_norm)
  
  # Add 'gene_name' column
  df_t$gene_id <- as.character(index_list[i])
  df_t$gene_name <- as.character(index_list_name[i])
  df_t$sample_id<- row.names(df_t)


  df_t <- merge(df_t, samplesheet[, c("sample_id", "condition", "cell_type")], by = "sample_id")
  # Filter for male
  baseline_df <- subset(df_t, condition == "10mix")
  
  # Filter for female
  perturbation_df <- subset(df_t, condition == "11mix")
  
  # Calculate annotation
  #anno = t.test(male_df$log1p_norm, female_df$log1p_norm)
  anno = wilcox.test(baseline_df$log1p_norm, perturbation_df$log1p_norm, alternative = "two.sided",  correct = TRUE)
  padjval = anno$p.value
  df_t$padjval = as.numeric(padjval)
  
  # Append the result
  result <- rbind(result, df_t)
  colnames(result)[2] <- "log1p_norm"
  
}


result <- result[!is.na(result$condition), ]
p  <-  ggviolin(result, x = "condition", y = "log1p_norm", trim=FALSE, title="Gene expression") +  
  facet_wrap(~ gene_name) 



max_y <- max(result$log1p_norm, na.rm = TRUE)

p <- p + geom_signif(comparisons = list(c("10mix", "11mix")),
                     map_signif_level = TRUE,
                     test = "t.test",
                     y_position = max_y + 1, 
                     tip_length = 0.02) # Adjust the y_position as needed
p
