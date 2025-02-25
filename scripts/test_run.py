import pandas as pd
import os
import requests
import time
from genes_extract import *


def submit_queries(cancer_name, method, num_genes, node_properties, network_properties, degs):
    genes_dict = extract_genes(
        method=method,
        num_genes=num_genes,
        node_properties_file=node_properties,
        network_properties_file=network_properties,
        degs_file=degs
    )

    upregulated_genes = genes_dict["upregulated"]
    downregulated_genes = genes_dict["downregulated"]

    if len(upregulated_genes) < 10 or len(downregulated_genes) < 10:
        st.error("Not enough genes extracted. Try increasing the gene count.")
        return None

    upregulated_genes_str = 'TAG\\t\\t' + '\\t'.join(upregulated_genes)
    downregulated_genes_str = 'TAG\\t\\t' + '\\t'.join(downregulated_genes)


    