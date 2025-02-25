library(biomaRt)
library(dplyr)
setwd("C:/Users/areeba khan/Documents/UdS/Master Thesis")
ABC_data <- read.csv("DE Analysis/protein_de_genes_DESEQ.csv")


ABC_data <- ABC_data %>% dplyr::filter(abs(log2FoldChange) > 1 )
hub_genes <- readLines("Tfmir Output/hub_genes.txt")
mds_genes <- readLines("mds_genes.txt")


########################### HUB GENES ###################
hub_genes <- as.data.frame(hub_genes)
colnames(hub_genes) <- c("gene")
filtered_hub_genes <- ABC_data[ABC_data$gene_name %in% hub_genes$gene, ]
head(filtered_hub_genes, 2)
# dim(filtered_hub_genes)

process_deseq2_data <- function(df, log2fc_threshold = 1) {
  # Filtering for up-regulated genes
  upregulated_genes <- df %>%
    dplyr::filter(log2FoldChange > log2fc_threshold) %>%
    arrange(pvalue)
  
  # Filtering for down-regulated genes
  downregulated_genes <- df %>%
    dplyr::filter(log2FoldChange < -log2fc_threshold) %>%
    arrange(pvalue)
  
  return(list(upregulated_genes = upregulated_genes, downregulated_genes = downregulated_genes))
}


filter <- function(df, n) {
  ## up regulated
  gene_ids_up <- df$upregulated_genes$gene_id
  ids_without_ver_up <- remove_version(gene_ids_up)
  entrez_upregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_up, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  upregulated_genes <- as.data.frame(entrez_upregulated$ENTREZID)
  top_n_up <- head(upregulated_genes, n)
  up_file_name <- paste0("Analysis/top_", n, "_hub_upregulated_genes.txt")
  write.table(top_n_up, up_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)

  ## down regulated
  gene_ids_down <- df$downregulated_genes$gene_id
  ids_without_ver_down <- remove_version(gene_ids_down)
  entrez_downregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_down, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  downregulated_genes <- as.data.frame(entrez_downregulated$ENTREZID)
  top_n_down <- head(downregulated_genes, n)
  down_file_name <- paste0("Analysis/top_", n, "_hub_downregulated_genes.txt")
  write.table(top_n_down, down_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)

}
results <- process_deseq2_data(filtered_hub_genes)
filter(results, 50)
filter(results, 75)
filter(results, 100)
filter(results, 150)

filter_and_separate_genes <- function(df, top_n) {
  top_genes <- df %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(top_n)
  
  #### up regulated
  upregulated_genes <- top_genes %>%
    dplyr::filter(log2FoldChange > 0)
  
  gene_ids_up <- upregulated_genes$gene_id
  ids_without_ver_up <- remove_version(gene_ids_up)
  entrez_upregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_up, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  upregulated_genes <- as.data.frame(entrez_upregulated$ENTREZID)
  up_file_name <- paste0("Analysis/variable_size/top_", top_n, "_upregulated_genes.txt")
  write.table(upregulated_genes, up_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
  

  ### down regulated
  downregulated_genes <- top_genes %>%
    dplyr::filter(log2FoldChange < 0)
  
  gene_ids_down <- downregulated_genes$gene_id
  ids_without_ver_down <- remove_version(gene_ids_down)
  entrez_downregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_down, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  downregulated_genes <- as.data.frame(entrez_downregulated$ENTREZID)
  down_file_name <- paste0("Analysis/variable_size/top_", top_n, "_downregulated_genes.txt")
  write.table(downregulated_genes, down_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)


}

filter_and_separate_genes(df, 100)
filter_and_separate_genes(df, 150)
filter_and_separate_genes(df, 200)
filter_and_separate_genes(df, 300)


####################################################### OLD FUNCTIONS ########################################################
top_n_equal <- function(df, n) {
  # Top N upregulated genes
  top_up <- head(df[order(df$pvalue, decreasing = TRUE), ], n)
  
  # Top N downregulated genes
  top_down <- head(df[order(df$pvalue, decreasing = FALSE), ], n)
  
#   return(rbind(top_up, top_down))
  return(list(top_up = top_up, top_down = top_down))
}

result <- top_n_equal(merge, 100)
top <- result$top_up
down <- result$top_down


dataframe_down <- merge %>% 
  filter(log2FoldChange <= -1) %>% 
  arrange(pvalue)

dataframe_up <- merge %>% 
  filter(log2FoldChange >= 1) %>% 
  arrange(pvalue)


top_100_down <- head(dataframe_down, 150)
top_100_down <- top_100_down[, c(2)]

top_100_up <- head(dataframe_up, 150)
top_100_up <- top_100_up[, c(2)]



write.table(top_100_down, "top_down_hub_150.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



############## MDS GENES ###################
mds_genes <- ens2symbol(mds_genes)
filtered_mds_genes <- ABC_data[ABC_data$gene_name %in% mds_genes$hgnc_symbol, ]
head(filtered_mds_genes, 2)

mds_genes_id <- ens2entrez(filtered_mds_genes$gene_name)
colnames(mds_genes_id) <- c("gene_name", "entrez_id")
mds_merge <- merge(mds_genes_id, filtered_mds_genes, by = "gene_name")


mds_result <- top_n_equal(mds_merge, 100)
mds_top <- mds_result$top_up
mds_down <- mds_result$top_down


dataframe_down_mds <- mds_merge %>% 
  filter(log2FoldChange <= -1) %>% 
  arrange(pvalue)

dataframe_up_mds <- mds_merge %>% 
  filter(log2FoldChange >= 1) %>% 
  arrange(pvalue)


top_100_down_mds <- head(dataframe_down_mds, 50)
top_100_down_mds <- top_100_down_mds[, c(2)]

top_100_up_mds <- head(dataframe_up_mds, 50)
top_100_up_mds <- top_100_up_mds[, c(2)]
df <- dataframe_down_mds[, c(2)]
write.table(df, "top_down_mds_all.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



################################################################################################

filtered_hub_genes$logFC_combined <- rowMeans(filtered_hub_genes[, c("logFC_A", "logFC_B", "logFC_C")], na.rm = TRUE)
ranked_hub_genes <- filtered_hub_genes[order(abs(filtered_hub_genes$logFC_combined), decreasing = TRUE), ]
top_hub_genes <- ranked_hub_genes[1:min(300, nrow(ranked_hub_genes)), ]
if (nrow(top_hub_genes) < 20) {
  top_hub_genes <- ranked_hub_genes[1:nrow(ranked_hub_genes), ]
} else {
  top_hub_genes <- ranked_hub_genes[1:300, ]
}

print(top_hub_genes)





top_genes <- head(ABC_data[order(ABC_data$logFC_A, decreasing = TRUE), ], 300) # Top 100 upregulated
bottom_genes <- head(ABC_data[order(ABC_data$logFC_A, decreasing = FALSE), ], 300) # Top 100 downregulated
selected_genes <- rbind(top_genes, bottom_genes)

selected_genes_id <- ens2entrez(bottom_genes$ensembl_gene_id)
list(selected_genes_id$entrezgene_id)
m <- merge(selected_genes_id, selected_genes, by = "ensembl_gene_id")
selected_genes_id


top_n_by_abs_logFC <- function(df, n) {
  return(head(df[order(abs(df$logFC_B), decreasing = TRUE), ], n))
}

#  DEGs by |log2(FC)|
top_150_degs <- top_n_by_abs_logFC(ABC_data, 200)

head(top_150_degs, 5)

write.csv(top_150_degs, "test.csv", row.names = FALSE)





top_n_equal <- function(df, n) {
  # Top N upregulated genes
  top_up <- head(df[order(df$logFC_C, decreasing = TRUE), ], n)
  
  # Top N downregulated genes
  top_down <- head(df[order(df$logFC_C, decreasing = FALSE), ], n)
  
  return(rbind(top_up, top_down))
}


top <- top_n_equal(ABC_data, 150)
g <- ens2entrez(top$ensembl_gene_id)
mer <- merge(g, top, by = "ensembl_gene_id")
write.csv(mer, "300_test.csv", row.names = FALSE)

protein_coding_de_genes <- read.csv("DE Analysis/mirna_deseq.csv")
significant_genes <- protein_coding_de_genes[abs(protein_coding_de_genes$log2FoldChange) > 1, ]

# Add a new column to represent regulation direction (1 for up, -1 for down)
significant_genes$regulation <- ifelse(significant_genes$log2FoldChange > 0, 1, -1)
significant_genes <- significant_genes[, c(1,3,6,7,8)]

de<-significant_genes[,c(1,5)]
write.table(de, file = "tfmir_mirna.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
