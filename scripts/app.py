import streamlit as st
import pandas as pd
import os
import requests
import time
import mygene
from genes_extract import *
from overlap_analysis import *
from io import BytesIO
import tarfile
import tempfile

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


st.title("CMAP Query Submission")
st.markdown("### Extract genes and submit queries to the CMAP API.")

cancer_name = st.text_input("Enter Cancer Type")
method = st.selectbox("Select Extraction Method", ["degs", "hub", "mds", "vertex", "bfs"])
num_genes = st.number_input("Number of Genes to Extract", min_value=10, max_value=150, value=100)

node_properties = st.file_uploader("Upload Node Properties File (TSV)", type=["tsv"], key="node_properties")
network_properties = st.file_uploader("Upload Network Properties File (YAML)", type=["yaml"], key="network_properties")
degs = st.file_uploader("Upload DEGs File", type=["csv"], key="degs")

if st.button("Submit Query"):
    if degs is None and (node_properties is None or network_properties is None):
        st.error("Please upload the DEGs file or both the node and network property files.")
    else:
        job_id = submit_queries(cancer_name, method, num_genes, node_properties, network_properties, degs)

        if job_id:
            st.session_state["job_id"] = job_id

def load_drug_data(file):
    if file is not None:
        return file.read().decode("utf-8").splitlines()
    else:
        st.info("No drug file uploaded. Using default drug file from Github.")
        try:
            response = requests.get()
            response.raise_for_status()
            return response.text.splitlines()
        except Exception as e:
            st.error(f"{e}")
            return []

def download_and_extract_gctx(job_id):
    headers = {
        'Accept': 'application/json',
        'user_key': '2f844a47b5bad670768300263f923a92',
    }

    # Step 1: Get job metadata
    info_res = requests.get(f'https://api.clue.io/api/jobs/findByJobId/{job_id}', headers=headers)
    if info_res.status_code != 200:
        st.error("Failed to fetch job metadata.")
        return None

    job_data = info_res.json()
    download_url = job_data.get('download_url', '')
    if download_url.startswith("//"):
        download_url = "https:" + download_url

    # Step 2: Download tar.gz file
    tar_res = requests.get(download_url)
    if tar_res.status_code != 200:
        st.error("Failed to download result archive.")
        return None

    tar_bytes = BytesIO(tar_res.content)

    # Step 3: Extract .gctx file to a temporary directory
    temp_extract_path = tempfile.mkdtemp()

    with tarfile.open(fileobj=tar_bytes, mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.name.endswith('/matrices/gutc/ps_pert_cell.gctx') and 'gutc/' in member.name:
                extracted_file = tar.extractfile(member)
                if extracted_file:
                    gctx_file_path = os.path.join(temp_extract_path, os.path.basename(member.name))
                    with open(gctx_file_path, 'wb') as out_file:
                        out_file.write(extracted_file.read())
                    return gctx_file_path

    st.warning("No GCTX file found in the archive.")
    return None


if st.button("Run Overlap"):
    drug_list = None 
    gctx_bytes = download_and_extract_gctx(st.session_state["job_id"])
    if gctx_bytes:
        cancer_name = cancer_name.strip()
        result = overlap_analysis_cmps(drug_list, gctx_bytes, cancer_name)
        st.write("Overlap Analysis Result:")
        st.dataframe(result)
    else:
        st.error("Failed to load GCTX data for overlap analysis.")


if "job_id" in st.session_state:
    if st.button("Fetch Results"):
        fetch_results(st.session_state["job_id"])