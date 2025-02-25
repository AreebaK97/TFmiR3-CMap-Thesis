import requests
import time
import csv
import os
from genes_extract import * 

def submit_queries(cancer_name, method, num_genes, node_properties_file=None, network_properties_file=None, degs_file=None):
    """
    Submit gene set queries to CMAP API.

    Parameters:
        cancer_name (str): Name of the cancer type.
        method (str): Method to extract genes ('degs', 'hub', 'mds', 'vertex', 'bfs').
        num_genes (int): Number of genes to extract.
        node_properties_file (str, optional): Path to the node properties file.
        network_properties_file (str, optional): Path to the network properties YAML file.
        degs_file (str, optional): Path to the DEGs file (required if method is 'degs').
    """
    
    dict_download = {}
    dict_errors = {}
    batch_size = 10
    query_count = 0
    query_name = f"{cancer_name}_{num_genes}_equal_{method}"

    try:
        # Extract genes using function instead of reading files
        genes_dict = extract_genes(
            method=method,
            num_genes=num_genes,
            node_properties_file=node_properties_file,
            network_properties_file=network_properties_file,
            degs_file=degs_file
        )

        upregulated_genes = genes_dict["upregulated"]
        downregulated_genes = genes_dict["downregulated"]

        if len(upregulated_genes) > 10 and len(downregulated_genes) > 10:
            upregulated_genes_str = 'TAG\t\t' + '\t'.join(upregulated_genes)
            downregulated_genes_str = 'TAG\t\t' + '\t'.join(downregulated_genes)

            headers = {
                'Content-Type': 'application/json',
                'Accept': 'application/json',
                'user_key': '2f844a47b5bad670768300263f923a92',
            }

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
                response_dict = response.json()
                dict_download[query_name] = response_dict['result']['job_id']
            else:
                dict_errors[query_name] = f"Error {response.status_code}: {response.text}"

            query_count += 1

            if query_count % batch_size == 0:
                print(f"Waiting for 30 minutes after submitting {query_count} queries...")
                time.sleep(1800)  # Pause for 30 minutes to prevent rate limits

    except Exception as e:
        dict_errors[query_name] = str(e)

    # Save results to CSV
    save_to_csv(dict_download, f'{cancer_name}_download_jobs.csv', ['query_name', 'job_id'])
    save_to_csv(dict_errors, f'{cancer_name}_error_jobs.csv', ['query_name', 'error'])


def save_to_csv(data_dict, file_name, fieldnames):
    """ Save dictionary data to a CSV file. """
    with open(file_name, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for key, value in data_dict.items():
            writer.writerow({fieldnames[0]: key, fieldnames[1]: value})


def download_files(job_ids_csv, output_directory, method_names):
    """
    Download query results from CMAP API.

    Parameters:
        job_ids_csv (str): Path to CSV file with job IDs.
        output_directory (str): Directory to save downloaded files.
        method_names (list): List of method names used.
    """

    headers = {
        'Accept': 'application/json',
        'user_key': '2f844a47b5bad670768300263f923a92',
    }

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    with open(job_ids_csv, mode='r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            query_name = row['query_name']
            job_id = row['job_id']
            response_id = requests.get(f'http://api.clue.io/api/jobs/findByJobId/{job_id}', headers=headers)

            if response_id.status_code == 200:
                di = response_id.json()
                download_url = di.get('download_url')

                if download_url and download_url.startswith("//"):
                    download_url = "https:" + download_url

                subfolder = determine_subfolder(query_name, method_names)
                save_path = os.path.join(output_directory, subfolder)

                if not os.path.exists(save_path):
                    os.makedirs(save_path)

                print(f"Fetching job result for: {query_name} (Job ID: {job_id})")

                response = requests.get(download_url, headers=headers)
                if response.status_code == 200:
                    local_filename = os.path.join(save_path, f"{query_name}.tar.gz")
                    print(f"Downloading {query_name} to {local_filename}...")
                    download_file(download_url, local_filename)
                    print(f"File saved as: {local_filename}")
                else:
                    print(f"Failed to fetch job details for Job ID: {job_id}. Status Code: {response.status_code}")

            else:
                print(f"Failed to retrieve job info for Job ID: {job_id}. Status Code: {response_id.status_code}")


def determine_subfolder(query_name, method_names):
    """
    Determine the subfolder to save files based on the query name.

    Parameters:
        query_name (str): Name of the query.
        method_names (list): List of method names.

    Returns:
        str: Path of the subfolder.
    """
    if "variable" in query_name.lower():
        for method in method_names:
            if method.lower() in query_name.lower():
                return os.path.join("variable_set", method)
        return "equal_set"
    
    for method in method_names:
        if method.lower() in query_name.lower():
            return method
    
    return "others"


def download_file(url, local_filename):
    """
    Download a file from a given URL and save it locally.

    Parameters:
        url (str): URL of the file.
        local_filename (str): Path to save the file.
    """
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(local_filename, 'wb') as file:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    file.write(chunk)
    except Exception as e:
        print(f"Error downloading file {local_filename}: {e}")



