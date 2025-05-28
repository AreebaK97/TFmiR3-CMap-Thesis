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
        degs_file: Path to the differentially expressed genes (DEGs) file (CSV format).
        log2fc_threshold (float, optional): Log2 fold change threshold for DEG filtering (default is 1).

    Returns:
        dict: A dictionary with 'upregulated' and 'downregulated' gene lists.
    """
    cancer_df = pd.read_csv(degs_file)
    output = {}
    if method == "degs":
        if not degs_file:
            raise ValueError("DEGs file must be provided when method is 'degs'.")

        pval_col = "padj" if "padj" in cancer_df.columns else "FDR"
        logfc_col = "log2FoldChange" if "log2FoldChange" in cancer_df.columns else "logFC"
        pval_raw_col = "pvalue" if "pvalue" in cancer_df.columns else "PValue"
        cancer_df = cancer_df[(cancer_df[pval_col] < 0.01) & (abs(cancer_df[logfc_col]) > log2fc_threshold)].sort_values(pval_raw_col)
        if "gene_name" in cancer_df.columns:
            genes = cancer_df["gene_name"].astype(str)
        else:
            genes = cancer_df.index.astype(str)

        genes = genes.str.replace(r"\.\d+$", "", regex=True)
        cancer_df = cancer_df.assign(cleaned_gene=genes)
        upregulated = cancer_df[cancer_df[logfc_col] > log2fc_threshold].head(num_genes)["cleaned_gene"].tolist()
        downregulated = cancer_df[cancer_df[logfc_col] < -log2fc_threshold].head(num_genes)["cleaned_gene"].tolist()
        
        # cancer_df = cancer_df[(cancer_df["padj"] < 0.01) & (abs(cancer_df["log2FoldChange"]) > log2fc_threshold)].sort_values("pvalue")
        #upregulated = cancer_df[cancer_df["log2FoldChange"] > log2fc_threshold].head(num_genes)["gene_name"].tolist()
        #downregulated = cancer_df[cancer_df["log2FoldChange"] < -log2fc_threshold].head(num_genes)["gene_name"].tolist()

    else:
        if not node_properties_file or not network_properties_file:
            raise ValueError("Both node properties and network properties files must be provided for this method.")

        node_data = pd.read_csv(node_properties_file, sep="\t")

        network_data = yaml.safe_load(network_properties_file)


        filtered_data = node_data[(node_data["node.type"] == "gene") | (node_data["node.type"] == "TF")]

        if method == "hub":
            genes = filtered_data.sort_values(by="degree", ascending=False)["node.ID"].tolist()
            method_genes = cancer_df[cancer_df["gene_name"].isin(genes)].sort_values(pval_raw_col)
            upregulated = method_genes[method_genes["log2FoldChange"] > log2fc_threshold].head(num_genes)["gene_name"].tolist()
            downregulated = method_genes[method_genes["log2FoldChange"] < -log2fc_threshold].head(num_genes)["gene_name"].tolist()
        elif method == "mds":
            genes = filtered_data[filtered_data["MDS"] == "True"]["node.ID"].tolist()
            method_genes = cancer_df[cancer_df["gene_name"].isin(genes)].sort_values(pval_raw_col)
            upregulated = method_genes[method_genes["log2FoldChange"] > log2fc_threshold].head(num_genes)["gene_name"].tolist()
            downregulated = method_genes[method_genes["log2FoldChange"] < -log2fc_threshold].head(num_genes)["gene_name"].tolist()
        elif method == "vertex":
            vertex_sort = network_data["hotspots"]["vertex.sort"]
            genes = list(vertex_sort["VS.core"] + vertex_sort["VS.top"])
            method_genes = cancer_df[cancer_df["gene_name"].isin(genes)].sort_values(pval_raw_col)
            upregulated = method_genes[method_genes["log2FoldChange"] > log2fc_threshold].head(num_genes)["gene_name"].tolist()
            downregulated = method_genes[method_genes["log2FoldChange"] < -log2fc_threshold].head(num_genes)["gene_name"].tolist()
        elif method == "bfs":
            bfs_dict = network_data["hotspots"]["BFS"]
            sorted_levels = sorted(bfs_dict.keys(), key=lambda x: int(x.split('.')[-1]))
            last_three_levels = sorted_levels[-3:]
            genes = []
            for level in last_three_levels:
                genes.extend(bfs_dict[level])
            method_genes = cancer_df[cancer_df["gene_name"].isin(genes)].sort_values(pval_raw_col)
            upregulated = method_genes[method_genes["log2FoldChange"] > log2fc_threshold].head(num_genes)["gene_name"].tolist()
            downregulated = method_genes[method_genes["log2FoldChange"] < -log2fc_threshold].head(num_genes)["gene_name"].tolist()
        else:
            raise ValueError("Invalid method! Choose from 'degs', 'hub', 'mds', 'vertex', 'bfs'.")

    output["upregulated"] = upregulated
    output["downregulated"] = downregulated
    return output

