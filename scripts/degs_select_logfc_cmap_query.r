library(data.table)
library(dplyr)
library(EnsDb.Hsapiens.v79)
setwd("C:/Users/areeba khan/Documents/UdS/Master Thesis/DE Analysis")

remove_version <- function(ids) {
  return(sub("\\..*$", "", ids))
}

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

df <- read.csv("protein_de_genes_DESEQ.csv")
df <- df %>% dplyr::filter(abs(log2FoldChange) > 1)

results <- process_deseq2_data(df)

remove_version <- function(ids) {
  return(sub("\\..*$", "", ids))
}

filter <- function(df, n) {
  ## up regulated
  gene_ids_up <- df$upregulated_genes$gene_id
  ids_without_ver_up <- remove_version(gene_ids_up)
  entrez_upregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_up, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  upregulated_genes <- as.data.frame(entrez_upregulated$ENTREZID)
  top_n_up <- head(upregulated_genes, n)
  up_file_name <- paste0("DEGS/equal_size/top_", n, "_upregulated_genes.txt")
  write.table(top_n_up, up_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)

  ## down regulated
  gene_ids_down <- df$downregulated_genes$gene_id
  ids_without_ver_down <- remove_version(gene_ids_down)
  entrez_downregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_down, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  downregulated_genes <- as.data.frame(entrez_downregulated$ENTREZID)
  top_n_down <- head(downregulated_genes, n)
  down_file_name <- paste0("DEGs/equal_size/top_", n, "_downregulated_genes.txt")
  write.table(top_n_down, down_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)

}


filter(results, 50)
filter(results, 75)
filter(results, 100)
filter(results, 150)


################################################ |log(FC)| Selected ######################################

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
  up_file_name <- paste0("DEGS/variable_size/top_", top_n, "_upregulated_genes.txt")
  write.table(upregulated_genes, up_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
  

  ### down regulated
  downregulated_genes <- top_genes %>%
    dplyr::filter(log2FoldChange < 0)
  
  gene_ids_down <- downregulated_genes$gene_id
  ids_without_ver_down <- remove_version(gene_ids_down)
  entrez_downregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_down, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  downregulated_genes <- as.data.frame(entrez_downregulated$ENTREZID)
  down_file_name <- paste0("DEGs/variable_size/top_", top_n, "_downregulated_genes.txt")
  write.table(downregulated_genes, down_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)


}

filter_and_separate_genes(df, 100)
filter_and_separate_genes(df, 150)
filter_and_separate_genes(df, 200)
filter_and_separate_genes(df, 300)
