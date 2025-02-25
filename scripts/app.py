import streamlit as st
import pandas as pd
import os
import requests
import time
import mygene
from genes_extract import *

def submit_queries(cancer_name, method, num_genes, node_properties, network_properties, degs):
    st.write("Extracting genes...")
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
    
    mg = mygene.MyGeneInfo()
    data = mg.querymany(upregulated_genes, scopes='symbol', fields='entrezgene', species='human')
    entrezgenelist = [item['entrezgene'] for item in data if 'entrezgene' in item]
    upregulated_genes_str = 'TAG\\t\\t' + ('\\t'.join([str(x) for x in entrezgenelist]))

    data = mg.querymany(downregulated_genes, scopes='symbol', fields='entrezgene', species='human')
    entrezgenelist = [item['entrezgene'] for item in data if 'entrezgene' in item]
    downregulated_genes_str = 'TAG\\t\\t' + ('\\t'.join([str(x) for x in entrezgenelist]))

    headers = {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'user_key': '2f844a47b5bad670768300263f923a92',
    }

    query_name = f"{cancer_name}_{num_genes}_{method}"
    data = f'''
    {{
        "tool_id": "sig_gutc_tool",
        "data_type": "L1000",
        "name": "{query_name}",
        "uptag-cmapfile": "{upregulated_genes_str}",
        "dntag-cmapfile": "{downregulated_genes_str}",
        "dataset": "Touchstone",
        "ignoreWarnings": true
    }}
    '''

    response = requests.post('https://api.clue.io/api/jobs', headers=headers, data=data)

    if response.status_code == 200:
        job_id = response.json()['result']['job_id']
        st.success(f"Query submitted successfully! Job ID: {job_id}")
        return job_id
    else:
        st.error(f"Error {response.status_code}: {response.text}")
        return None


# Function to fetch and download results
def fetch_results(job_id):
    headers = {
        'Accept': 'application/json',
        'user_key': '2f844a47b5bad670768300263f923a92',
    }
    response = requests.get(f'http://api.clue.io/api/jobs/findByJobId/{job_id}', headers=headers)

    if response.status_code == 200:
        job_data = response.json()
        download_url = job_data.get('download_url', '')

        if download_url.startswith("//"):
            download_url = "https:" + download_url

        st.success("Job completed! Click below to download results.")
        st.markdown(f"[Download Results]({download_url})")

    else:
        st.error(f"Error fetching job results. Status Code: {response.status_code}")


# Streamlit UI
st.title("CMAP Query Submission")
st.markdown("### Extract genes and submit queries to the CMAP API.")

# User Inputs
cancer_name = st.text_input("Enter Cancer Type", "prad")
method = st.selectbox("Select Extraction Method", ["degs", "hub", "mds", "vertex", "bfs"])
num_genes = st.number_input("Number of Genes to Extract", min_value=10, max_value=500, value=100)

node_properties = st.file_uploader("Upload Node Properties File (TSV)", type=["tsv"], key="node_properties")
network_properties = st.file_uploader("Upload Network Properties File (YAML)", type=["yaml"], key="network_properties")
degs = st.file_uploader("Upload DEGs File (if using DEGs method)", type=["csv"], key="degs")

if st.button("Submit Query"):
    if degs is None and (node_properties is None or network_properties is None):
        st.error("Please upload at least the DEGs file or both the node and network property files.")
    else:
        job_id = submit_queries(cancer_name, method, num_genes, node_properties, network_properties, degs)

        if job_id:
            st.session_state["job_id"] = job_id

if "job_id" in st.session_state:
    if st.button("Fetch Results"):
        fetch_results(st.session_state["job_id"])