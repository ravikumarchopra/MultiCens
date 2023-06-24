import os
import re
import json
import gzip
import time
import asyncio
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
import statsmodels.api as sm
from pyensembl import EnsemblRelease
from statsmodels.tools.tools import add_constant
from sklearn.linear_model import LinearRegression

# Load the Ensembl database
ensembl_db = EnsemblRelease(86)

# Dataset paths
gtex_path = "dataset/GTEx_Analysis_v8_eQTL_expression_matrices/"
gtex_cov_path = "dataset/GTEx_Analysis_v8_eQTL_covariates/"
other_path = "upload/"
other_cov_path = "upload/"


# Map GENEID to SYMBOL
def get_gene_names(ensembl_gene_ids):
    # Map GENEID to SYMBOL
    gene_mapping = []
    for ensembl_gene_id in ensembl_gene_ids:
        try:
            gene = ensembl_db.gene_by_id(ensembl_gene_id)
            if gene:
                gene_mapping.append(
                    {"GENEID": ensembl_gene_id, "SYMBOL": gene.gene_name}
                )
        except Exception as e:
            # print(str(e))
            pass

    return pd.DataFrame(gene_mapping)

# Remove rows which does not have a SYMBOL
def remove_extra_ensembl_ids(tissue_df, HGNC_gene_mapping):
    # Remove rows which does not have a SYMBOL
    merged_df = tissue_df.merge(
        HGNC_gene_mapping, left_on="gene_id", right_on="GENEID", how="left"
    )
    merged_df = merged_df[merged_df["SYMBOL"].notnull()]
    symbols = merged_df["SYMBOL"]
    merged_df.drop(columns=["GENEID", "SYMBOL"], inplace=True)
    return merged_df, symbols

def add_column(tissue_incomplete, column_name):
    num_genes = tissue_incomplete.shape[0]
    random_vector = np.random.uniform(low=0.0001, high=0.00011, size=num_genes)
    df = pd.DataFrame(tissue_incomplete)
    df[column_name] = random_vector
    tissue_complete = df.rename(columns={"random_vector": column_name})
    return tissue_complete.to_numpy()

# Read gene expression data for a particular tssiue from GTEx data
def get_tissue(tissue):
    # Read gene expression data for a tissue from GTEx data
    try:
        if tissue["dataset"] == "GTex":
            tissue_df = pd.read_csv(tissue["exp"], sep="\t", compression="gzip")
        else:
            tissue_df = pd.read_csv(tissue["exp"])
    except FileNotFoundError:
        print("Tissue file not found")
        return None
    tissue_df["gene_id"] = tissue_df["gene_id"].str[:15]
    ensembl_gene_list = tissue_df["gene_id"]
    # Map GENEID to SYMBOL
    HGNC_gene_mapping = get_gene_names(ensembl_gene_list)
    # Remove extra ensembl ids for which there is no SYMBOL
    tissue_df, headers = remove_extra_ensembl_ids(tissue_df, HGNC_gene_mapping)
    gene_tpm = tissue_df.iloc[:, 5:].T
    gene_tpm.columns = headers
    return gene_tpm

# Remove the effect of co-variates
def rm_cov_effect(tissue_name, gene_tpm, cov_path):
    try:
        # Reading covariates for the tissue
        covariates = pd.read_csv(cov_path, sep="\t").set_index("ID")
        covariates = covariates.T
    except FileNotFoundError:
        print("Covariates file not found")
        return None

    headers = gene_tpm.columns
    # gene_tpm.columns = [col.replace(".", "-") for col in gene_tpm.columns]
    common_samples = gene_tpm.index.intersection(covariates.index)
    gene_tpm = gene_tpm.loc[common_samples].T
    covariates = covariates.loc[common_samples]

    residuals = []
    print(
        f"\n------------------------------ [ {tissue_name} ] ------------------------------\n"
    )
    # Remove the effect of co-variates
    for i, row in tqdm(
        gene_tpm.iterrows(),
        total=len(gene_tpm),
        desc="Genes Processed: ",
        unit="Gene(s)",
    ):
        # Perform linear regression to remove covariate effects
        model = sm.OLS(row, sm.add_constant(covariates))
        result = model.fit()
        # Get the residuals
        residuals.append(result.resid)

    tissue = pd.DataFrame(residuals).T
    tissue.columns = headers

    return tissue

# Create correlation matrix
def create_corr_matrix(tissue_list, session_id):
    start_time = time.time()
    n = len(tissue_list)
    tissue_names = "_".join(sorted([t["name"] for t in tissue_list]))
    output_path = os.path.join("uploads", session_id)
    corr_file_path = os.path.join(output_path, f"{tissue_names}.csv")
    if os.path.exists(corr_file_path):
        print("\n[Correlation matrix for selected tissues already exists.]")
        print("\n[Loading correlation matrix...]")
        correlation_matrix = pd.read_csv(corr_file_path)
        correlation_matrix = correlation_matrix.set_index(correlation_matrix.columns[0])
        with open(os.path.join(output_path, f"{tissue_names}.json"), 'r') as file:
            gene_count = json.loads(file.read())
    else:  
        tissues, sample_ids, genes = [], [], []
        # Read gene expression data for each tissue
        print("\n[Reading gene expression data for each tissue]")
        print("-------------------------------------------------")
        for i in tqdm(range(n), total=n, desc="Tissues Processed: ", unit=" Tissue(s)"):
            tissue = get_tissue(tissue_list[i])
            tissues.append(tissue)
            sample_ids.append(set(tissue.index))
            genes.append(set(tissue.columns))

        common_samples = list(set.intersection(*sample_ids))
        common_genes = list(set.intersection(*genes))
        print(f'\nCommon Samples: {len(common_samples)}')
        print(f'\nCommon Genes: {len(common_genes)}')
        most_var_genes = set()
        # Find 2000 most varying genes for each tissue
        print("\n[Find 2000 most varying genes for each tissue]")
        print("-------------------------------------------------")
        for i in tqdm(range(n), total=n, desc="Tissues Processed: ", unit=" Tissue(s)"):
            # Select Common samples
            tissues[i] = tissues[i].loc[common_samples]
            tissues[i] = tissues[i].loc[:, ~tissues[i].columns.duplicated()]
            # Select common genes
            tissues[i] = tissues[i][common_genes]
            # Calculate variances for tissue columns
            variances = np.var(tissues[i], axis=0)
            # Sort source_variances in descending order and get the top 2000 indices
            sorted_indices = np.argsort(variances)[::-1][:2000]
            sorted_columns = tissues[i].columns[sorted_indices]
            most_var_genes.update(sorted_columns)
            # print(len(sorted_columns))
        
        most_var_genes = list(most_var_genes)
        gene_count = []
        # Time calculation
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time/60, 2)
        time_taken = f'{execution_time} min(s)'
        print("\nTime taken to find the common samples and genes and select top 2K genes from each tissue:", time_taken, "\n")
        
        # Remove the effect of co-variates
        print("\n[Removing the effect of co-variates]")
        print("-------------------------------------------------")
        start_time = time.time()
        for i in range(n):
            # Save gene list for a Tissue
            # genes = pd.DataFrame(tissues[i].columns)
            # genes.to_csv(
            #     os.path.join(
            #         "static", "downloads", tissue_list[i]["name"].capitalize() + ".csv"
            #     ),
            #     header=False,
            #     index=False,
            # )
            # Select the union of (5000 most varying genes from each tissue) for each tissue
            tissues[i] = tissues[i][most_var_genes]
            gene_count.append(len(tissues[i].columns))
            if tissue_list[i]["dataset"] == "GTex":
                tissues[i] = rm_cov_effect(
                    tissue_list[i]["name"].capitalize(), tissues[i], tissue_list[i]["cov"]
                )
            tissues[i].columns = [f"{gene}.{i+1}" for gene in tissues[i].columns]
        
        # Time calculation
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time/60, 2)
        time_taken = f'{execution_time} min(s)'
        print("\nTime taken to remove the effect of covariates from each tissue:", time_taken, "\n")        
        
        start_time = time.time()
        # Concatenating tissues data horizontally
        tissues_df = pd.concat(tissues, axis=1, ignore_index=False)

        # Calculate correlation matrix
        print("\n[Creating correlation matrix]")
        print("-------------------------------------------------")
        with tqdm(total=len(tissues_df.columns), desc="Column-pairs Processed: ", unit=" Column-pair(s)") as pbar:
            def update_pbar(x):
                pbar.update()
                return x
            correlation_matrix = tissues_df.corr(method="spearman", min_periods=1).apply(
                update_pbar, axis=1
            )
        # Convert correlation score to absolute numbers
        correlation_matrix = correlation_matrix.abs()
        
        # Time calculation
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time/60, 2)
        time_taken = f'{execution_time} min(s)'
        print("\nTime taken to concatenate the tissues horizontally and computing the correlation matrix:", time_taken, "\n")
        
        # Saving correlation matrix to a csv file
        os.makedirs(output_path, exist_ok=True)
        # asyncio.run(save_correlation_matrix_async(correlation_matrix, gene_count, output_path, tissue_names))
        thread = threading.Thread(target=save_correlation_matrix_async, args=(correlation_matrix, gene_count, output_path, tissue_names))
        thread.start()

    return correlation_matrix, gene_count

def save_correlation_matrix_async(correlation_matrix, gene_count, output_path, tissue_names):
    correlation_matrix.to_csv(
            os.path.join(output_path, f"{tissue_names}.csv")
        )
    with open(os.path.join(output_path, f"{tissue_names}.json"), 'w') as file:
            json.dump(gene_count, file)
    print('[Correlation matrix saved.]')

def allowed_file(filename):
    ALLOWED_EXTENSIONS = {"zip", "gz"}
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS
