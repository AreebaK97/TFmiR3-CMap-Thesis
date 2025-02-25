library(dplyr)
library(EnsDb.Hsapiens.v79)
library(utils)
library(yaml)
setwd("C:/Users/areeba khan/Documents/UdS/Master Thesis")

# ABC_data <- read.csv("DE Analysis/protein_de_genes_DESEQ.csv")
# ABC_data <- ABC_data %>% dplyr::filter(abs(log2FoldChange) > 1 ) # significant genes selected only those that have abs(logFC) greater than 1
# hub_genes <- readLines("Tfmir Output/hub_genes.txt")
# mds_genes <- readLines("Tfmir Output/mds_genes.txt")
# vertex_genes <- readLines("Tfmir Output/vertex_sort_genes.txt")
# bfs_genes <- readLines("Tfmir Output/bfs_genes.txt")

cancers <- c("brca", "luad", "lusc", "prad", "lihc", "coad", "skcm")
methods <- c("degs", "hub", "mds", "vertex", "bfs")


create_method_files <- function(cancer_tfmir_path, cancer_name, file){
  
  node_properties_path <- paste0(cancer_tfmir_path, "/results/network-complete/node_properties.tsv")
  node_properties_path_yml <- paste0(cancer_tfmir_path, "/results/network-complete/network_properties.yaml")

  node_data <- read.delim(node_properties_path, sep = "\t")
  node_data_yml <- yaml.load_file(node_properties_path_yml)

  filtered_data <- dplyr::filter(node_data, node.type == 'gene' | node.type == 'TF')

  hub_genes <- filtered_data %>%
    arrange(desc(degree)) %>%
    pull(node.ID)
  

  mds_genes <- dplyr::filter(filtered_data, MDS == "True") %>% 
    pull(node.ID)



  vertex_sort <- node_data_yml$hotspots$'vertex.sort'
  vertex_core <- vertex_sort$'VS.core'
  vertex_top <- vertex_sort$'VS.top'


  bfs_genes <- node_data_yml$hotspots$BFS
  bfs_level0 <- bfs_genes$'BFS.level.0'
  bfs_level1 <- bfs_genes$'BFS.level.1'
  

  hub_output <- paste0(file, "/hub_genes.txt")
  mds_output <- paste0(file, "/mds_genes.txt")
  vertex_output <- paste0(file, "/vertex_genes.txt")
  bfs_output <- paste0(file, "/bfs_genes.txt")


  writeLines(hub_genes, con = hub_output)
  writeLines(mds_genes, con = mds_output)

  writeLines(c(
    Filter(function(gene) !grepl('hsa', gene), vertex_core),
    Filter(function(gene) !grepl('hsa', gene), vertex_top)),con = vertex_output)

  writeLines(c(
    Filter(function(gene) !grepl('hsa', gene), bfs_level0),
    Filter(function(gene) !grepl('hsa', gene), bfs_level1) ), con = bfs_output)
}

for (cancer in cancers){
  tryCatch({file_path <- paste0(getwd(),"/tfmir_output/", cancer, "_tfmir_output.zip")
  unzip_file_path <- paste0(getwd(),"/tfmir_output/", cancer, "_tfmir_output")
  if (!file.exists(file_path)) {
      stop("Zip file does not exist:", file_path)
    }
  unzip(file_path, exdir = unzip_file_path)
  maindir <- paste0(getwd(), "/gene_sets/selected_genes/", cancer)
  if (dir.exists(maindir)){
    create_method_files(unzip_file_path,cancer, maindir)
  } else {
    dir.create(maindir, recursive = TRUE)
    create_method_files(unzip_file_path, cancer, maindir)
  }}, error = function(e){
    cat("Failed to Proccess ", toupper(cancer), " Data: ", conditionMessage(e), "\n\n")
  })
  
}

gene_num <- list(100, 150, 200, 300)

for (cancer in cancers){
  for (method in methods){

    tryCatch({
      cancer_degs <- process_deseq2_data(cancer, method)
      significant_genes <- cancer_degs$genes
      cancer_genes <- cancer_degs$df
    }, error = function(e) {
      cat("Failed to Proccess ", toupper(cancer), " Data: ", conditionMessage(e), "\n\n")
      next
    })

    for (n in gene_num){
      tryCatch({
        num = n / 2
        filter(significant_genes, num, cancer, method)
      }, error = function(e){
        cat("Failed to Filter Equal Gene Sets For ", toupper(cancer), " Data, Method: ",toupper(method), " For ", num, " Gene Set. Error: ", conditionMessage(e), "\n\n")
      })


      tryCatch({
        filter_and_separate_genes(cancer_genes, n, cancer, method)
      }, error = function(e){
        cat("Failed to Filter Variable Gene Sets For ", toupper(cancer), " Data, Method: ", toupper(method), " For ", n, " Gene Set. Error: ", conditionMessage(e), "\n\n")
      })
    }
  }
}


remove_version <- function(ids) {
  return(sub("\\..*$", "", ids))
}

process_deseq2_data <- function(cancer_name, method, log2fc_threshold = 1) {
  # Read Cancer DE Analysis File
  cancer_deg_file <- paste0("DE Analysis/de_results/", cancer_name, "_de_genes_deseq2_padj.csv")
  cancer_df <- read.csv(cancer_deg_file)
  cancer_df <- cancer_df %>% dplyr::filter(padj < 0.01) # significant genes selection based on adjusted p-value
  cancer_df <- cancer_df %>% dplyr::filter(abs(log2FoldChange) > 1 ) # significant genes selected only those that have abs(logFC) greater than 1

  if (method == "degs"){
    df <- cancer_df
  } else{
    # Read Method Selected (Hubs, MDS, BFS, Vertex Sort) Files for Each Cancer
    method_selected_genes_file <- paste0("gene_sets/selected_genes/", cancer_name, "/", method, "_genes.txt")
    method_genes <- readLines(method_selected_genes_file)
    method_genes <- as.data.frame(method_genes)
    colnames(method_genes) <- c("gene")
    df <- cancer_df[cancer_df$gene_name %in% method_genes$gene, ]
  }


  # Filtering for up-regulated genes
  upregulated_genes <- df %>%
    dplyr::filter(log2FoldChange > log2fc_threshold) %>%
    arrange(pvalue)
  
  # Filtering for down-regulated genes
  downregulated_genes <- df %>%
    dplyr::filter(log2FoldChange < -log2fc_threshold) %>%
    arrange(pvalue)
  
  return(list(df = df, genes = list(upregulated_genes = upregulated_genes, downregulated_genes = downregulated_genes)))
}


#### EQUAL GENE SET
filter <- function(df, n, cancer_name, method) {

  # checking the existence of the directory
  maindir <- paste0(getwd(), "/gene_sets/cmap_query_sets/", cancer_name)
  if (dir.exists(maindir)){
    subdir <- paste0(maindir, "/equal_sets")
    dir.create(subdir, recursive = TRUE)
    subdir <- paste0(subdir, "/", method)
    dir.create(subdir, recursive = TRUE)
  } else {
    dir.create(maindir, recursive = TRUE)
    subdir <- paste0(maindir, "/equal_sets")
    dir.create(subdir, recursive = TRUE)
    subdir <- paste0(subdir, "/", method)
    dir.create(subdir, recursive = TRUE)
  }

  ## up regulated
  gene_ids_up <- df$upregulated_genes$gene_id
  ids_without_ver_up <- remove_version(gene_ids_up)
  entrez_upregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_up, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  upregulated_genes <- as.data.frame(entrez_upregulated$ENTREZID)
  top_n_up <- head(upregulated_genes, n)
  up_file_name <- paste0(subdir, "/top_", n, "_upregulated_genes.txt")
  write.table(top_n_up, up_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)

  ## down regulated
  gene_ids_down <- df$downregulated_genes$gene_id
  ids_without_ver_down <- remove_version(gene_ids_down)
  entrez_downregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_down, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  downregulated_genes <- as.data.frame(entrez_downregulated$ENTREZID)
  top_n_down <- head(downregulated_genes, n)
  down_file_name <- paste0(subdir, "/top_", n, "_downregulated_genes.txt")
  write.table(top_n_down, down_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)

}


#### VARIABLE GENE SET
filter_and_separate_genes <- function(df, top_n, cancer_name, method) {

  # checking the existence of the directory
  maindir <- paste0(getwd(), "/gene_sets/cmap_query_sets/", cancer_name)
  if (dir.exists(maindir)){
    subdir <- paste0(maindir, "/variable_sets")
    dir.create(subdir, recursive = TRUE)
    subdir <- paste0(subdir, "/", method)
    dir.create(subdir, recursive = TRUE)
  } else {
    dir.create(maindir, recursive = TRUE)
    subdir <- paste0(maindir, "/variable_sets")
    dir.create(subdir, recursive = TRUE)
    subdir <- paste0(subdir, "/", method)
    dir.create(subdir, recursive = TRUE)
  }
  
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
  up_file_name <- paste0(subdir, "/top_", top_n, "_upregulated_genes.txt")
  write.table(upregulated_genes, up_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
  

  ### down regulated
  downregulated_genes <- top_genes %>%
    dplyr::filter(log2FoldChange < 0)
  
  gene_ids_down <- downregulated_genes$gene_id
  ids_without_ver_down <- remove_version(gene_ids_down)
  entrez_downregulated <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids_without_ver_down, keytype = "GENEID", columns = c("SYMBOL","GENEID", "ENTREZID"))
  downregulated_genes <- as.data.frame(entrez_downregulated$ENTREZID)
  down_file_name <- paste0(subdir, "/top_", top_n, "_downregulated_genes.txt")
  write.table(downregulated_genes, down_file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}







#####################################################################
########################### HUB GENES ###############################
hub_genes <- as.data.frame(hub_genes)
colnames(hub_genes) <- c("gene")
filtered_hub_genes <- ABC_data[ABC_data$gene_name %in% hub_genes$gene, ]
results <- process_deseq2_data(filtered_hub_genes)

### Equal gene set size
filter(results, 50, 'DE Analysis/Hub/equal_size')
filter(results, 75, 'DE Analysis/Hub/equal_size')
filter(results, 100, 'DE Analysis/Hub/equal_size')
filter(results, 150, 'DE Analysis/Hub/equal_size')

### Variable gene set size
filter_and_separate_genes(filtered_hub_genes, 100, "DE Analysis/Hub/variable_size")
filter_and_separate_genes(filtered_hub_genes, 150, "DE Analysis/Hub/variable_size")
filter_and_separate_genes(filtered_hub_genes, 200, "DE Analysis/Hub/variable_size")
filter_and_separate_genes(filtered_hub_genes, 300, "DE Analysis/Hub/variable_size")


#####################################################################
########################### MDS GENES ###############################
mds_genes <- as.data.frame(mds_genes)
colnames(mds_genes) <- c("gene")
filtered_mds_genes <- ABC_data[ABC_data$gene_name %in% mds_genes$gene, ]
results_mds <- process_deseq2_data(filtered_mds_genes)

### Equal gene set size
filter(results_mds, 25, 'DE Analysis/MDS/equal_size')
filter(results_mds, 50, 'DE Analysis/MDS/equal_size')
# filter(results_mds, 100, 'DE Analysis/MDS/equal_size')
# filter(results_mds, 150, 'DE Analysis/MDS/equal_size')

### Variable gene set size
filter_and_separate_genes(filtered_mds_genes, 50, "DE Analysis/MDS/variable_size")
filter_and_separate_genes(filtered_mds_genes, 100, "DE Analysis/MDS/variable_size")
filter_and_separate_genes(filtered_mds_genes, 150, "DE Analysis/MDS/variable_size")




#######################################################################
########################### Vertex Sort ###############################
vertex_genes <- as.data.frame(vertex_genes)
colnames(vertex_genes) <- c("gene")
filtered_vertex_genes <- ABC_data[ABC_data$gene_name %in% vertex_genes$gene, ]
results_vertex <- process_deseq2_data(filtered_vertex_genes)

### Equal gene set size
filter(results_vertex, 50, "DE Analysis/Vertex_Sort/equal_size")
filter(results_vertex, 75, "DE Analysis/Vertex_Sort/equal_size")
filter(results_vertex, 100, "DE Analysis/Vertex_Sort/equal_size")
filter(results_vertex, 150, "DE Analysis/Vertex_Sort/equal_size")
### Variable gene set size
filter_and_separate_genes(filtered_vertex_genes, 100, "DE Analysis/Vertex_Sort/variable_size")
filter_and_separate_genes(filtered_vertex_genes, 150, "DE Analysis/Vertex_Sort/variable_size")
filter_and_separate_genes(filtered_vertex_genes, 200, "DE Analysis/Vertex_Sort/variable_size")
filter_and_separate_genes(filtered_vertex_genes, 300, "DE Analysis/Vertex_Sort/variable_size")



########################################################################
########################### BFS SELECTED ###############################
bfs_genes <- as.data.frame(bfs_genes)
colnames(bfs_genes) <- c("gene")
filtered_bfs_genes <- ABC_data[ABC_data$gene_name %in% bfs_genes$gene, ]
results_bfs <- process_deseq2_data(filtered_bfs_genes)

### Equal gene set size
filter(results_bfs, 50, "DE Analysis/BFS/equal_size")
filter(results_bfs, 75, "DE Analysis/BFS/equal_size")
filter(results_bfs, 100, "DE Analysis/BFS/equal_size")
filter(results_bfs, 150, "DE Analysis/BFS/equal_size")
### Variable gene set size
filter_and_separate_genes(filtered_bfs_genes, 100, "DE Analysis/BFS/variable_size")
filter_and_separate_genes(filtered_bfs_genes, 150, "DE Analysis/BFS/variable_size")
filter_and_separate_genes(filtered_bfs_genes, 200, "DE Analysis/BFS/variable_size")
filter_and_separate_genes(filtered_bfs_genes, 300, "DE Analysis/BFS/variable_size")



######################################################
############### Volcano Plot #########################



library(readr)
library(ggplot2)

df <- read_csv("DE Analysis/mirna_deseq.csv")

alpha <- 0.01

df$`-log10(padj)` <- -log10(df$padj)

df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]

ggplot(df, aes(x = log2FoldChange, y = `-log10(padj)`)) +
  geom_point(aes(color = padj < 0.01), alpha = 0.3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme_minimal()
