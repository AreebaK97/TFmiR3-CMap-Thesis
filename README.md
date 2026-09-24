# TFmiR3 and CMap: A Workflow for Drug Prioritization
This repository contains the computational workflow and Streamlit application developed as part of my Master's thesis, TFmiR3 and CMap: A Workflow for Drug Prioritization, at Universität des Saarlandes.

The project integrates differential expression analysis, mRNA-miRNA co-regulatory networks (TFmiR3), network-based gene prioritization and the Connectivity Map (CMap) to identify and evaluate potential cancer drug candidates.

A Python-based Streamlit App was developed to automate the downstream drug discovery workflow and provide an interactive interface for submitting refined gene signatures to the CMap Query App API. 
The app can be accessed freely at https://tfmir3-cmap.streamlit.app/

## Project Objective

The objective of this workflow was to develop and validate a computational pipeline for drug repurposing in cancer treatment by integrating gene regulatory network analysis (using TFmiR3: freely available at https://service.bioinformatik.uni-saarland.de/tfmir3-test/) with the Connectivity Map (CMap) database. The goal was to improve drug discovery efficiency by using transcriptional and network-based methods to identify and prioritize existing compounds with potential therapeutic value for cancer treatment.

## Workflow

TCGA Data

  &darr;

Differential Expression Analysis

  &darr;

TFmiR3 Co-regulatory Network Construction

  &darr;

Network-based Gene Prioritization

  &darr;

Refined Gene Signatures

  &darr;

CMap Drug Query

  &darr;

Candidate Drug Validation  

## Streamlit App

A Python-based Streamlit App was developed to automate the downstream CMap drug discovery workflow.

### Inputs
The application requires:
- Differential expression analysis results (.csv)
- TFmiR3 output files i.e. network_properties.yaml and node_properties.tsv
- Cancer type being analysed
- Network analysis method (Hub, MDS, Vertex Sort or BFS)
- CMap API key from CLUE account

If no network analysis method is selected, the application uses the differentially expressed genes as the input signature.

For drug validation, users can either use the available FDA-approved list for the selected cancer or provide a custom drug list.

### Application Workflow
**1. Signature Extraction and Formatting:** Extract key genes from TFmiR3 outputs based on the selected network method; otherwise, use differently expressed genes.
Convert gene symbols to ENTREZ IDs, as the CMap API only works with BING Entrez IDs; other identifiers will not work.

**2. Query Submission Automation:** Formats the input into GMT (Gene Matrix Transposed) files, which are then submitted to the CMap server via a POST request using the user’s API key. Once analysis is completed, the app uses the job ID to send a request to the CMap API to download the results, which is a compressed TAR file containing GTCX matrix files. The file ps_pert_cell.gctx stores compounds and the associated cell-line-specific connectivity scores.

**3. Drug Validation:** Compares predicted compounds with FDA-approved drugs for the selected cancer type, or with a user-provided custom drug list (using PubChempy for compound ID mapping). Performs overlap analysis to assess the reliability of the predictions.

**4. Result Output:** Compares predicted compounds with FDA-approved drugs for the selected cancer type, or with a user-provided custom drug list (using PubChempy for compound ID mapping).Performs overlap analysis to assess the reliability of the predictions.
