import tarfile
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from cmapPy.pandasGEXpress.parse import parse
from functools import reduce
import json
import pubchempy as pcp
import requests
from io import StringIO

cancer_celllines = {"brca" : "MCF7", "luad": "A549", "lusc" : "A549", "lihc" : "HEPG2", "prad" : "PC3", "coadread" : "HT29"}

def read_gctx_file(file, cancer_name):
    cmp_data = parse(file)
    metadata = cmp_data.row_metadata_df
    cmap_scores = cmp_data.data_df
    # cmapdrugs_scores
    cell_line = cancer_celllines[cancer_name]
    metadata.reset_index(drop = True, inplace = True)
    mcf7_brd_df = metadata[metadata['pert_id'].str.contains('BRD') & metadata['cell_id'].str.contains(cell_line) & metadata['pert_type'].str.contains('trt_cp')]
    
    cmap_scores.reset_index(inplace=True)
    cmap_scores.columns = ["cid", "TAG"]

    mcf7_scores = cmap_scores[cmap_scores['cid'].str.contains('BRD') & cmap_scores['cid'].str.contains(cell_line)]
    # mcf7_scores['pert_id'] = mcf7_scores['cid'].str.split(':', expand=True)[0]
    mcf7_scores.loc[:, 'pert_id'] = mcf7_scores['cid'].str.split(':', expand=True)[0]

    merged_df = pd.merge(mcf7_scores, mcf7_brd_df[['pert_id', 'pert_iname', 'pert_type']], on='pert_id')
    final_df = merged_df[['pert_id', 'pert_iname', 'TAG']]
    cmap_name_cids = pd.read_csv(f"https://raw.githubusercontent.com/AreebaK97/TFmiR3-CMap-Thesis/main/cids/{cancer_name}/{cancer_name}_cmap_cmps_cid.txt", header=None, sep="\t")

    cmap_name_cids.columns = ['pert_iname', 'cid']

    cmap_name_cids['pert_iname'] = cmap_name_cids['pert_iname'].str.lower()

    cmap_name_cids.drop_duplicates(inplace=True)

    cmap_name_cids = cmap_name_cids.dropna(subset=['cid'])

    cmap_name_cids = cmap_name_cids.groupby('pert_iname')['cid'].agg(lambda x: ','.join(x.astype(str).unique())).reset_index()
    cmap_name_cids['pert_iname'] = cmap_name_cids['pert_iname'].str.lower()

    final_df['pert_iname'] = final_df['pert_iname'].str.lower()
    cmap = pd.merge(cmap_name_cids, final_df[['pert_id', 'pert_iname', 'TAG']], on='pert_iname')

    cmap['pert_iname'] = cmap['pert_iname'].str.lower()

    cmap_grouped = cmap.groupby(['pert_iname', 'TAG', 'pert_id'], as_index=False).agg({
        'cid': lambda x: ','.join(x)
    })

    cmap = cmap_grouped.drop_duplicates()

    def unique_join(val):
        cids = val.split(',')
        cids = [str(int(float(cid))) for cid in cids]
        return ','.join(sorted(set(cids)))

    cmap['cid'] = cmap['cid'].apply(unique_join)

    cmap.reset_index(drop=True, inplace=True)

    return cmap

def drug_filter(cid_file, ref_drug):
    response = requests.get(cid_file)
    response.raise_for_status()

    file_like = StringIO(response.text)

    lines = file_like.readlines()
    drugs = []
    ids = []

    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            drugs.append(parts[0].lower())
            ids.append(parts[1])

    df = pd.DataFrame({"Drug": drugs, "Cid": ids})

    # df.dropna(inplace=True)
    # df.reset_index(drop=True, inplace=True)
    df = df.groupby('Drug')['Cid'].agg(lambda x: ','.join(x.astype(str).unique())).reset_index()
    df['Drug'] = df['Drug'].str.lower()
    df_merge = pd.merge(df, ref_drug["Drug"], on = "Drug")
    return df_merge


def remove_drug(drug_df, cmap):
    keep = []
    for i in range(len(drug_df)):
        cids = drug_df.iloc[i, 1].split(",") 
        for k in range(len(cmap)):
            cids_cmap = str(cmap.iloc[k, 3])
            # print(cids_cmap)
            if isinstance(cids_cmap, str):
                cids_cmap = cids_cmap.split(",")
                if any(cid in cids_cmap for cid in cids):
                    keep.append(i)
                    break

    filtered_drug = drug_df.iloc[keep]
    return filtered_drug

def get_compound_ids(drug_name):
    try:
        compounds = pcp.get_compounds(drug_name, 'name')
        
        compound_ids = [compound.cid for compound in compounds]
        
        return compound_ids
    except Exception as e:
        print("An error occurred:", e)
        return []
    
def function_call(drug_list, file, cancer_name):    
    cmap = read_gctx_file(file, cancer_name)
    drug_cid = {}
    # print(cmap)
    if drug_list is None:
        ttd_drug = pd.read_csv(f"https://raw.githubusercontent.com/AreebaK97/TFmiR3-CMap-Thesis/main/Drugs/{cancer_name}/ttd_drugs.txt", 
            sep="\t", header=None, names=["Drug", "ID"])
        ttd_drug.iloc[:, 0] = ttd_drug.iloc[:, 0].str.lower()
        ttd_drug.drop_duplicates(inplace=True)

        nci_ref_drugs = pd.read_csv(f"https://raw.githubusercontent.com/AreebaK97/TFmiR3-CMap-Thesis/main/Drugs/{cancer_name}/nci_drugs.txt", 
            sep="\t", header=None, names=["Drug"])
        nci_ref_drugs.iloc[:, 0] = nci_ref_drugs.iloc[:, 0].str.lower()
        nci_ref_drugs.drop_duplicates(inplace=True)

        ttd_drug_cids = drug_filter(f"https://raw.githubusercontent.com/AreebaK97/TFmiR3-CMap-Thesis/main/Drugs/{cancer_name}/ttd_cid.txt", 
            ttd_drug)
        nci_drug_cids = drug_filter(f"https://raw.githubusercontent.com/AreebaK97/TFmiR3-CMap-Thesis/main/Drugs/{cancer_name}/nci_cid.txt", 
            nci_ref_drugs)


        nci = remove_drug(nci_drug_cids, cmap)
        ttd = remove_drug(ttd_drug_cids, cmap)

        approved_drugs = pd.concat([ttd, nci], ignore_index=True)
        approved_drugs['Drug'] = approved_drugs['Drug'].str.lower()
        approved_drugs = approved_drugs.groupby('Drug')['Cid'].apply(lambda x: ','.join(x)).reset_index()

        approved_drugs.drop_duplicates(inplace=True)

        approved_drugs['Cid'] = approved_drugs['Cid'].apply(lambda x: ','.join(set(x.split(','))))
    
    else:
        for drug in drug_list:
            if "(" in drug: 
                parts = drug.split('(')
                drug_name1 = parts[0].strip()  # Name outside parentheses
                drug_name2 = parts[1].split(')')[0].strip()  # Name inside parentheses

                compound_ids1 = get_compound_ids(drug_name1)
                compound_ids2 = get_compound_ids(drug_name2)

                if compound_ids1 or compound_ids2:
                    if compound_ids1 in compound_ids2 or compound_ids2 in compound_ids1:
                        if len(compound_ids1) == len(compound_ids2):
                            if drug in drug_cid.keys():
                                drug_cid[drug] = compound_ids1
                            else:
                                drug_cid[drug] = compound_ids1
                    else:
                        if len(compound_ids1) > len(compound_ids2):
                            if drug_name1 in drug_cid.keys():
                                drug_cid[drug_name1] = compound_ids1
                            else:
                                drug_cid[drug_name1] = compound_ids1
                        if len(compound_ids1) > len(compound_ids2):
                            if drug_name2 in drug_cid.keys():
                                drug_cid[drug_name2] = compound_ids2
                            else:
                                drug_cid[drug_name2] = compound_ids2
            else:
                compound_ids = get_compound_ids(drug)
                # print("Compound IDs for", drug, ":", compound_ids)
                if compound_ids:
                    if drug in drug_cid.keys():
                        drug_cid[drug] = compound_ids
                    else:
                        drug_cid[drug] = compound_ids
    keep = []
    for i in range(len(approved_drugs)):
        cids = set(approved_drugs.iloc[i, 1].split(','))
        for k in range(len(approved_drugs)):
            cids_cmap = set(approved_drugs.iloc[k, 1].split(','))
            cut = cids.intersection(cids_cmap)
            if cut and len(cids) < len(cids_cmap):
                keep.append(i)
                break

    if keep:
        approved_drugs = approved_drugs.drop(keep).reset_index(drop=True)
    
    cmap_grouped_sorted = cmap.sort_values(by='TAG', ascending=True).reset_index(drop=True)
    #print(approved_drugs)
    
    return cmap_grouped_sorted, approved_drugs


def overlap_coefficient(compound_df, approved_drugs_df, n):
    count = 0
    first_n_items = compound_df.head(n)
    # print(first_n_items)
    # print(approved_drugs_df)
    overlapping_compounds = []    
    for i, row in first_n_items.iterrows():
        cids = str(row['cid']).split(",")
        cmp_name = row["pert_iname"]
        pert_id = row["pert_id"]
        # print(cids)
        for j, approved_row in approved_drugs_df.iterrows():
            cids_cmap = str(approved_row['Cid']).split(",")
            drug_name = approved_row["Drug"]
            # if any(cid in cids_cmap for cid in cids):
            #     count += 1
            #     break
            overlap = [cid for cid in cids if cid in cids_cmap]
            if overlap:
                count += 1
            
                overlapping_compounds.append({
                    "cmap_compunds": cmp_name,
                    "drug" : drug_name,
                    "overlap_cid" : overlap,
                    "pert_id" : pert_id
                })
    
    min_size = min(len(first_n_items), len(approved_drugs_df))
    # print(count)
    coeff = count / min_size
    # print(overlapping_compounds)
    return round(coeff * 100), overlapping_compounds


def overlap_coefficient_relevant(compound_df, approved_drugs_df):
    # Filter compounds with TAG value <= -95
    filtered_compound_df = compound_df[compound_df['TAG'] <= -95]
    
    count = 0
    
    for i, row in filtered_compound_df.iterrows():
        cids = str(row['cid']).split(",")
        for j, approved_row in approved_drugs_df.iterrows():
            cids_cmap = str(approved_row['Cid']).split(",")
            if any(cid in cids_cmap for cid in cids):
                count += 1
                break
    
    min_size = min(len(filtered_compound_df), len(approved_drugs_df))
    if min_size > 0:
        coeff = count / min_size
    else:
        coeff = 0
    return round(coeff * 100)


def overlap_analysis_cmps(drug_list, gtcx_file_path, cancer):

    cmap_mds, approved_drugs = function_call(drug_list, gtcx_file_path, cancer)
    # print(cmap_mds)
    relevant_score = overlap_coefficient_relevant(cmap_mds, approved_drugs)

    unique_overlap_cmps = set()
    all_overlap_cmps = []
    for n in [10, 20, 50, 100, 150, 200]:
        score, overlaped_cmps = overlap_coefficient(cmap_mds, approved_drugs, n)
        # print(overlaped_cmps)
        all_overlap_cmps.extend(overlaped_cmps)
    return(overlaped_cmps)
