import pandas as pd
import yaml


def extract_genes(method, num_genes, node_properties_file=None, network_properties_file=None, degs_file=None, log2fc_threshold=1):
    """
    Extract genes based on the specified method and number of genes.

    Parameters:
        method (str): The method to use ('degs', 'hub', 'mds', 'vertex', 'bfs').
        num_genes (int): Number of genes to extract.
        node_properties_file (str, optional): Path to the node properties file (TSV format).
        network_properties_file (str, optional): Path to the network properties file (YAML format).
        degs_file (str, optional): Path to the differentially expressed genes (DEGs) file (CSV format).
        log2fc_threshold (float, optional): Log2 fold change threshold for DEG filtering (default is 1).

    Returns:
        dict: A dictionary with 'upregulated' and 'downregulated' gene lists.
    """

    output = {}
    if method == "degs":
        if not degs_file:
            raise ValueError("DEGs file must be provided when method is 'degs'.")

        cancer_df = pd.read_csv(degs_file)
        cancer_df = cancer_df[(cancer_df["padj"] < 0.01) & (abs(cancer_df["log2FoldChange"]) > log2fc_threshold)].sort_values("pvalue")

        # upregulated = cancer_df[cancer_df["log2FoldChange"] > log2fc_threshold].nlargest(num_genes, "pvalue")["gene_name"].tolist()
        # downregulated = cancer_df[cancer_df["log2FoldChange"] < -1].nsmallest(num_genes, "pvalue")["gene_name"].tolist()
        upregulated = cancer_df[cancer_df["log2FoldChange"] > log2fc_threshold].head(num_genes)["gene_name"].tolist()

        downregulated = cancer_df[cancer_df["log2FoldChange"] < -log2fc_threshold].head(num_genes)["gene_name"].tolist()

    else:
        if not node_properties_file or not network_properties_file:
            raise ValueError("Both node properties and network properties files must be provided for this method.")

        node_data = pd.read_csv(node_properties_file, sep="\t")

        with open(network_properties_file, "r") as file:
            network_data = yaml.safe_load(file)

        filtered_data = node_data[(node_data["node.type"] == "gene") | (node_data["node.type"] == "TF")]

        if method == "hub":
            genes = filtered_data.sort_values(by="degree", ascending=False)["node.ID"].head(num_genes).tolist()
        elif method == "mds":
            genes = filtered_data[filtered_data["MDS"] == "True"]["node.ID"].head(num_genes).tolist()
        elif method == "vertex":
            vertex_sort = network_data["hotspots"]["vertex.sort"]
            genes = (vertex_sort["VS.core"] + vertex_sort["VS.top"])[:num_genes]
        elif method == "bfs":
            bfs_genes = network_data["hotspots"]["BFS"]
            genes = (bfs_genes["BFS.level.0"] + bfs_genes["BFS.level.1"])[:num_genes]
        else:
            raise ValueError("Invalid method! Choose from 'degs', 'hub', 'mds', 'vertex', 'bfs'.")

        upregulated, downregulated = genes[:num_genes], []

    output["upregulated"] = upregulated
    output["downregulated"] = downregulated
    return output

# genes_dict = extract_genes(method="degs", num_genes=100, degs_file="C:/Users/areeba khan/Documents/UdS/Master Thesis/DE Analysis/de_results/brca_de_genes_deseq2_padj.csv")


# upregulated_genes = genes_dict["upregulated"]
# downregulated_genes = genes_dict["downregulated"]

# if len(upregulated_genes) < 10 or len(downregulated_genes) < 10:
#     print("TOO SMALL")


# import mygene

# mg = mygene.MyGeneInfo()
# data = mg.querymany(upregulated_genes, scopes='symbol', fields='entrezgene', species='human')
# entrezgenelist = [item['entrezgene'] for item in data if 'entrezgene' in item]

# upregulated_genes_str = 'TAG\\t\\t' + ('\\t'.join([str(x) for x in entrezgenelist]))
# downregulated_genes_str = 'TAG\\t\\t' + '\\t'.join(downregulated_genes)

# print(upregulated_genes_str)

# print("\n\n")

# cancer_names = ["brca"]
# methods = ["degs"]
# n = 100
# for cancer in cancer_names:
#         for method in methods:
#                 new_n = n + n
#                 query_name = cancer + "_" + str(new_n) + "_equal_" + method
#                 with open(f"C:/Users/areeba khan/Documents/UdS/Master Thesis/gene_sets/cmap_query_sets/{cancer}/equal_sets/{method}/top_{n}_upregulated_genes.txt") as file:
#                     file_list = file.read().splitlines()

#                 with open(f"C:/Users/areeba khan/Documents/UdS/Master Thesis/gene_sets/cmap_query_sets/{cancer}/equal_sets/{method}/top_{n}_downregulated_genes.txt") as f:
#                     down_f = f.read().splitlines()

#                 if len(down_f) > 10 and len(file_list) > 10:
#                     file_list = 'TAG\\t\\t' + ('\\t'.join([str(x) for x in file_list]))
#                     down_f = 'TAG\\t\\t' + ('\\t'.join([str(x) for x in down_f]))

                

#                     print((file_list))

# # submit_queries("brca", "degs", 100, degs="C:/Users/areeba khan/Documents/UdS/Master Thesis/DE Analysis/de_results/brca_de_genes_deseq2_padj.csv")

