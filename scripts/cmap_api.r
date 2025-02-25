library(cluequery)
setwd("C:/Users/areeba khan/Documents/UdS/Master Thesis/gene_sets/cmap_query_sets/")

CLUE_API_KEY = "2f844a47b5bad670768300263f923a92"

cancers <- c("prad", "coad")
methods <- c("degs", "hub", "mds", "vertex", "bfs")

gene_set_num <- list(150, 200, 300)

for (cancer in cancers){
  for (method in methods){
    queries <- list()
    for (n in gene_set_num){
      up_genes_file_path <- paste0(cancer, "/variable_sets/", method, "/top_", n, "_upregulated_genes.txt")
      down_genes_file_path <- paste0(cancer, "/variable_sets/", method, "/top_", n, "_downregulated_genes.txt")

      if (file.exists(up_genes_file_path) && file.exists(down_genes_file_path)){
        up_gene <- readLines(up_genes_file_path)
        down_gene <- readLines(down_genes_file_path)
        name_cm <- paste0(n, "_", cancer, "_", method, "_variable_set")
        pre_gmt <- clue_gmt_from_list(up_gene, down_gene, name = name_cm)
        queries[[name_cm]] <- list(up = pre_gmt[["up"]], down = pre_gmt[["down"]])
        # query_submit <- clue_query_submit(pre_gmt[["up"]], pre_gmt[["down"]], name = name_cm, api_key = CLUE_API_KEY)
        # clue_query_wait(query_submit, interval = 60, timeout = 5000, api_key = CLUE_API_KEY)

      } else {
        cat("Skipping", cancer, method, "with set size", n, "as one or both files are missing.\n\n")
      }
    }
    queries_submit <- clue_queries_submit(queries = queries, api_key = CLUE_API_KEY, interval = 60)
    # clue_query_wait(queries_submit, interval = 60, timeout = 5000, api_key = CLUE_API_KEY)

    for (query_name in names(queries_submit)) {
      job_id <- queries_submit[[query_name]]$result$job_id
      clue_query_wait(job_id, interval = 60, timeout = 5000, api_key = CLUE_API_KEY)
    }

    for (query_name in names(queries_submit)) {
      parts <- unlist(strsplit(query_name, "_"))
      gene_set <- parts[1]  
      cancer <- parts[2]
      method <- parts[3]    

      new_query_name <- paste(cancer, gene_set, method, sep = "_")
      job_id <- queries_submit[[query_name]]$result$job_id
      destination_path <- paste0("C:/Users/areeba khan/Documents/UdS/Master Thesis/CMAP_Output/", cancer, "/variable_set/", method, "/", new_query_name, ".tar.gz")      
      clue_query_download(clue_query = job_id, destination = destination_path, api_key = CLUE_API_KEY)
      cat("Downloaded results for", query_name, "to", destination_path, "\n") 
    }
  }
}

for (query_name in names(queries_submit)) {
    job_id <- queries_submit[[query_name]]$result$job_id
    destination_path <- paste0("C:/Users/areeba khan/Documents/UdS/Master Thesis/CMAP_Output/", cancer, "/equal_set/hub/", query_name, ".tar.gz")      
    clue_query_download(clue_query = job_id, destination = destination_path, api_key = CLUE_API_KEY)
    cat("Downloaded results for", query_name, "to", destination_path, "\n") }


c <- clue_list_jobs(api_key = CLUE_API_KEY)
c$created
#####################################################################################################################

down <- readLines("prad/variable_sets/degs/top_100_downregulated_genes.txt")
up <- readLines("prad/variable_sets/degs/top_100_upregulated_genes.txt")

pre_gmt <- clue_gmt_from_list(up, down, "test")
sb <- clue_query_submit(pre_gmt[["up"]], pre_gmt[["down"]], name = "testing_api", api_key = CLUE_API_KEY)
clue_query_wait(sb, api_key = "2f844a47b5bad670768300263f923a92")
outputfile <- "mds_query.tar.gz"
result_path <- clue_query_download(query_submit, destination = outputfile, api_key = "2f844a47b5bad670768300263f923a92")


lusc_vertex_200 <- readLines("lusc/variable_sets/vertex/top_200_upregulated_genes.txt")







lusc_vertex_200_d <- readLines("lusc/variable_sets/vertex/top_200_downregulated_genes.txt")
pre_gmt <- clue_gmt_from_list(lusc_vertex_200, lusc_vertex_200_d, "test")




name <- "lusc_vertex"
up_gmt <- paste(name, "up", paste(lusc_vertex_200, collapse = "\t"), sep = "\t")
down_gmt <- paste(name, "down", paste(lusc_vertex_200_d, collapse = "\t"), sep = "\t")


writeLines(up_gmt, "upprad.gmt")
writeLines(down_gmt, "downprad.gmt")


require ( httr )
query <-  'example_query'
post_url <- 'https://api.clue.io/api/jobs'

files <- list(
    'tool_id' = 'sig_gutc_tool',
    'data_type' = 'L1000',
    "dataset" = "Touchstone",
    'name' = name,
    'up-cmapfile' = upload_file("upprad.gmt"),
    'dntag-cmapfile' = upload_file("downprad.gmt"),
    'ignoreWarnings' = TRUE
  )
  
headers <- c('user_key' = CLUE_API_KEY)

res <- POST(
  url = post_url,
  add_headers(.headers = headers),
  body = files)

contents <- content(res)
jID <- contents$result$job_id








prad_100_vertex <- clue_gmt_from_list(prad_vertex_100, prad_vertex_100_d, name = "100_prad_bfs")
prad_150_vertex <- clue_gmt_from_list(prad_vertex_150, prad_vertex_150_d, name = "150_prad_bfs")
prad_200_vertex <- clue_gmt_from_list(prad_vertex_200, prad_vertex_200_d, name = "200_prad_bfs")
prad_300_vertex <- clue_gmt_from_list(prad_vertex_300, prad_vertex_300_d, name = "300_prad_bfs")



queries <- list(
  prad_bfs_100 = list(
    up = prad_100_vertex[["up"]],
    down = prad_100_vertex[["down"]]
  ),
  prad_bfs_150 = list(
    up = prad_150_vertex[["up"]],
    down = prad_150_vertex[["down"]]
  ),
  prad_bfs_200 = list(
    up = prad_200_vertex[["up"]],
    down = prad_200_vertex[["down"]]
  ),
  prad_bfs_300 = list(
    up = prad_300_vertex[["up"]],
    down = prad_300_vertex[["down"]]
  )
)

ds <- clue_queries_submit(queries = queries, api_key = CLUE_API_KEY, interval = 60)
job_id <- ds$prad_bfs_100$result$job_id

outputfile <- "vertex_query.tar.gz"
result_path <- clue_query_download(job_id, destination = outputfile, api_key = CLUE_API_KEY)

clue_query_poll(ds, api_key = CLUE_API_KEY)
r <- clue_parse_result(
   result_path,
   score_level =  "cell",
   score_type = "tau",
   result_type ="pert"
 )

head(r$pert_id) 
typeof(r)



library(dplyr)
library(stringr)
mcf <- r %>% filter(str_detect(id, "BRD") & str_detect(id, "MCF7"))
mcf
dim(mcf)
mcf_fil <- mcf %>% filter(tau <= -95)
mcf_fil %>% arrange(tau)

clue_retrieve_api_key()



library(httr)
library(jsonlite)
url <- "https://api.clue.io/api/jobs"
headers <- add_headers(
  `Content-Type` = "application/json",
  `Accept` = "application/json",
  `user_key` = "2f844a47b5bad670768300263f923a92"
)

body <- list(
  "tool_id" = "sig_gutc_tool",
  "data_type" = "L1000",
  "name" = "Pitavastatin treated HUVEC cells (1 uM at 4H) vs. DMSO treated",
  "uptag-cmapfile" = upload_file("upprad.gmt"),
  "dntag-cmapfile" = upload_file("downprad.gmt"),
  "dataset" = "Touchstone"
  )

response <- POST(url, headers, body = body)

content(response, "text")

















library(httr)

# Load gene lists from text files
up_genes <- readLines("lusc/variable_sets/vertex/top_200_upregulated_genes.txt")
down_genes <- readLines("lusc/variable_sets/vertex/top_200_downregulated_genes.txt")

uptag_cmapfile <- paste("Upregulated_Genes", "", paste(up_genes, collapse = "\t"), sep = "\t")
dntag_cmapfile <- paste("Downregulated_Genes", "", paste(down_genes, collapse = "\t"), sep = "\t")

data <- list(
  tool_id = "sig_gutc_tool",
  data_type = "L1000",
  name = "(GSE32547) Pitavastatin treated HUVEC cells (1 uM at 4H) vs. DMSO treated",
  uptag_cmapfile = "TAG\t\t10365\t1831\t9314\t4846\t678\t22992\t3397\t26136\t79637\t5551\t7056\t79888\t1032\t51278\t64866\t29775\t994\t51696\t81839\t23580\t219654\t57178\t7014\t57513\t51599\t55818\t4005\t4130\t4851\t2050\t50650\t9469\t54438\t3628\t54922\t3691\t65981\t54820\t2261\t2591\t7133\t162427\t10912\t8581\t2523\t25807\t9922\t30850\t4862\t8567\t79686\t55615\t51283\t3337\t2887\t3223\t6915\t6907\t26056\t259217\t6574\t23097\t5164\t57493\t7071\t5450\t113146\t8650",
  dntag_cmapfile = "TAG\t\t5128\t5046\t956\t10426\t9188\t23403\t7204\t1827\t3491\t9076\t330\t8540\t22800\t10687\t19\t63875\t10979\t51154\t10370\t50628\t7128\t6617\t7187\t22916\t81034\t58516\t3096\t4794\t5202\t26511\t8767\t2355\t22943\t1490\t133\t11010\t51025\t23160\t56902\t3981\t5209\t6347\t5806\t7357\t9425\t3399\t6446\t64328\t6722\t8545\t688\t861\t390\t23034\t51330\t51474\t2633\t4609",
  dataset = "Touchstone"
)

headers <- add_headers(
  `Content-Type` = "application/json",
  `Accept` = "application/json",
  user_key = CLUE_API_KEY  # Replace with your actual API key
)

# Make the POST request
response <- POST("https://api.clue.io/api/jobs", body = toJSON(data), encode = "json", headers)

# Check the response
print(response$status_code)
print(content(response, "text"))

