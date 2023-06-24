import os
import re
import time
import uuid
import json
import shutil
import itertools
import numpy as np
import pandas as pd
from pathlib import Path
from memory_profiler import profile
from util.ranking import getRanking
from werkzeug.utils import secure_filename
from datetime import datetime, timedelta, timezone
from util.util import allowed_file, create_corr_matrix
from flask import Flask, render_template, request, session, redirect, flash
from util.centrality import (
    global_centrality,
    local_centrality,
    right_target_global_centrality_t,
)

tissue_names = {
    "pituitary": "Pituitary",
    "adrenal": "Adrenal_Gland",
    "adipose": "Adipose_Visceral_Omentum",
    "muscle": "Muscle_Skeletal",
    "heart": "Heart_Atrial_Appendage",
    "kidney": "Kidney_Cortex",
    "liver": "Liver",
    "thyroid": "Thyroid",
    "hypothalamus": "Brain_Hypothalamus",
    "ovary": "Ovary",
    "intestine": "Small_Intestine_Terminal_Ileum",
    "stomach": "Stomach",
    "pancreas": "Pancreas",
    "uterus": "Uterus",
    "breast": "Breast_Mammary_Tissue",
}

tissue_list = [
    "adipose",
    "liver",
    "muscle",
    "intestine",
    "thyroid",
    "ovary",
    "breast",
    "hypothalamus",
    "kidney",
    "pancreas",
    "pituitary",
    "adrenal",
    "uterus",
    "stomach",
    "heart",
]
hormone_list = [
    "adrenaline",
    "aldosterone",
    "angiotensin",
    "cortisol",
    "estradiol",
    "glucagon",
    "insulin",
    "leptin",
    "norepinephrine",
    "progesterone",
    "somatotrophin",
    "thyroxin",
    "vitamind",
]
measures = {
    "local": "Local Centrality",
    "global": "Global Centrality",
    "query": "Query Set Centrality",
}

dataset_list = ["GTex", "Others"]
gtex_path = "dataset/GTEx_Analysis_v8_eQTL_expression_matrices"
gtex_cov_path = "dataset/GTEx_Analysis_v8_eQTL_covariates"

tissue_list.sort()
hormone_list.sort()
# tissue_list.append('other')

app = Flask(__name__)
app.secret_key = "$#ghFkGH56-1aRGFGHGtrfJH"
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_TYPE"] = "filesystem"
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(
    minutes=60
)  # Session timeout of 1 hour
session_id = ""


# Session management
@app.before_request
def before_request():
    try:
        check_session_expiry()
    except Exception:
        pass


def check_session_expiry():
    if "last_activity" in session:
        last_activity_time = session["last_activity"]
        session_timeout = app.config["PERMANENT_SESSION_LIFETIME"]
        # Check if session has expired
        if datetime.now(timezone.utc) - last_activity_time > session_timeout:
            shutil.rmtree(os.path.join("uploads", session_id))
            file_path = os.path.join("downloads", session_id)
            if os.path.exists(file_path):
                os.remove(file_path)
            session.clear()

    # Update last activity time
    session["last_activity"] = datetime.now(timezone.utc)


# Home page
@app.route("/")
def index():
    global session_id
    session_id = str(uuid.uuid4().hex)
    session["id"] = session_id
    return render_template(
        "home.html",
        tissue_list=tissue_list,
        hormone_list=hormone_list,
        measures=measures,
        dataset_list=dataset_list,
    )


# Tool v1.0
@app.route("/get_ranking", methods=["POST"])
def get_ranking():
    if request.method == "POST":
        start_time = time.time()
        if len(request.files) <= 0:
            flash("No selected file")
            return redirect(request.url)

        file1 = request.files["input1"]
        file2 = request.files["input2"]
        source_tissue = request.form.get("source_tissue")
        target_tissue = request.form.get("target_tissue")
        hormone = request.form.get("hormone")

        source_path = "static/uploads/" + uuid.uuid4().hex + Path(file1.filename).suffix
        target_path = "static/uploads/" + uuid.uuid4().hex + Path(file2.filename).suffix
        file1.save(source_path)
        file2.save(target_path)

        # result = pd.read_csv(source)
        result = pd.DataFrame()

        try:
            result = getRanking(
                source_path, target_path, hormone, source_tissue, target_tissue
            )
        except:
            pass

        # delete file(s) from the server
        os.remove(source_path)
        os.remove(target_path)

        result["rank"] = range(1, len(result.index) + 1, 1)
        result.columns = ["Gene Name", "Centrality", "Rank"]

        result_csv = "static/downloads/" + uuid.uuid4().hex + ".csv"

        result.to_csv(result_csv, index=False)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time / 60, 2)
        # unit = 'sec(s)' if execution_time<1 else 'min(s)'
        # time_taken = f'{execution_time} {unit}'
        time_taken = f"{execution_time} min(s)"
        print("Execution time:", time_taken, "\n")

        return render_template(
            "ranking.html",
            columns=result.columns,
            rows=list(result.values.tolist()),
            path=result_csv,
            time_taken=time_taken,
        )


@app.route("/run_sample", methods=["GET"])
def run_sample():
    if request.method == "GET":
        start_time = time.time()
        source_tissue = "pancreas"
        target_tissue = "muscle"
        hormone = "insulin"

        source_path = "static/downloads/insulin_source_genes_full.csv"
        target_path = "static/downloads/insulin_target_genes_full.csv"

        result = pd.DataFrame()

        try:
            result = getRanking(
                source_path, target_path, hormone, source_tissue, target_tissue
            )
        except Exception as e:
            print(str(e))

        result["rank"] = range(1, len(result.index) + 1, 1)
        result.columns = ["Gene Name", "Centrality", "Rank"]

        result_csv = "static/downloads/" + uuid.uuid4().hex + ".csv"

        result.to_csv(result_csv, index=False)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time / 60, 2)
        # unit = 'sec(s)' if execution_time<1 else 'min(s)'
        # time_taken = f'{execution_time} {unit}'
        time_taken = f"{execution_time} min(s)"
        print("Execution time:", time_taken, "\n")

        return render_template(
            "ranking.html",
            columns=result.columns,
            rows=list(result.values.tolist()),
            path=result_csv,
            time_taken=time_taken,
        )


# Tool v2.0 (beta)
@app.route("/example_run", methods=["POST"])
@profile
def example_run():
    if request.method == "POST":
        start_time = time.time()
        measure = request.form.get("measure")
        tissues = json.loads(request.form.get("tissues"))
        # tissues = sorted(tissues, key=lambda t: t['name'])
        result = pd.DataFrame()

        try:
            for tissue in tissues:
                tissue_name = tissue_names[tissue["name"].lower()]
                tissue["exp"] = os.path.join(
                    gtex_path, f"{tissue_name}.v8.normalized_expression.bed.gz"
                )
                tissue["cov"] = os.path.join(
                    gtex_cov_path, f"{tissue_name}.v8.covariates.txt"
                )

            A, gene_count = create_corr_matrix(tissues, session_id)
            n = len(A) / len(tissues)
            if measure == "query":
                query_tissue = request.form.get("query-tissue")
                query_file = os.path.join("static", "downloads", "Pancreas.csv")
                query_index = int(request.form.get("tissue-index"))
                n = len(A) / len(tissues)
                s = int(query_index * n)
                e = int(s + n)
                query_gene_list = pd.read_csv(query_file, header=None)[0].tolist()
                query_gene_list = [
                    f"{gene}.{query_index+1}" for gene in query_gene_list
                ]
                genes = A.columns[s:e]
                common_target_genes = np.intersect1d(genes, query_gene_list)
                genes_indices = [
                    i for i, e in enumerate(A.columns) if e in common_target_genes
                ]
                _, g = right_target_global_centrality_t(
                    A.values,
                    num_layers=len(tissues),
                    target_tissue=query_index,
                    target_gene_indices=genes_indices,
                    start=s,
                    end=e,
                    p=0.9,
                )
                result["Centrality"] = g
            elif measure == "global":
                result["Centrality"] = global_centrality(A.values, len(tissues), p=0.9)
            else:
                result["Centrality"] = local_centrality(A.values, len(tissues), p=0.9)
            name_list = [
                [tissues[i]["name"] for _ in range(gene_count[i])]
                for i in range(len(tissues))
            ]
            result["Tissue Name"] = list(itertools.chain(*name_list))
            result["Gene Name"] = [re.sub(r"\.\d+$", "", s) for s in A.columns]
            if measure == "local":
                ranking = pd.DataFrame(
                    columns=["Tissue Name", "Gene Name", "Centrality", "Rank"]
                )
                for tissue in tissues:
                    df = result.loc[result["Tissue Name"] == tissue["name"]]
                    df = df.sort_values(by="Centrality", ascending=False)
                    df["Rank"] = [i + 1 for i in range(len(df))]
                    ranking = pd.concat([ranking, df], ignore_index=True)
                result = ranking
            else:
                result = result.sort_values(by="Centrality", ascending=False)
                result["Rank"] = [i + 1 for i in range(len(result))]
            result = result[["Tissue Name", "Gene Name", "Centrality", "Rank"]]
        except Exception as e:
            print(str(e))

        result = result.loc[result["Centrality"] != 0]
        result_csv = os.path.join("static", "downloads", session_id + ".csv")
        result.to_csv(result_csv)
        result["Centrality"] = np.round(result["Centrality"], 4)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time / 60, 2)
        # unit = 'sec(s)' if execution_time<1 else 'min(s)'
        time_taken = f"{execution_time} min(s)"
        print("Execution time:", time_taken, "\n")

        return render_template(
            "ranking.html",
            columns=result.columns,
            rows=list(result.values.tolist()),
            path=result_csv,
            time_taken=time_taken,
        )


@app.route("/get_centrality", methods=["POST"])
@profile
def get_centrality():
    if request.method == "POST":
        start_time = time.time()
        measure = request.form.get("measure")
        tissues = json.loads(request.form.get("tissues"))
        # tissues = sorted(tissues, key=lambda t: t['name'])
        result = pd.DataFrame()
        try:
            folder_path = os.path.join("uploads", session_id)
            os.makedirs(folder_path, exist_ok=True)
            files = request.files.getlist("files")
            for file in files:
                if file and allowed_file(file.filename):
                    filename = secure_filename(file.filename)
                    file_path = os.path.join(folder_path, filename)
                    file.save(file_path)
            for tissue in tissues:
                tissue_name = tissue_names[tissue["name"].lower()]
                if tissue["dataset"] == "GTex":
                    tissue["exp"] = os.path.join(
                        gtex_path, f"{tissue_name}.v8.normalized_expression.bed.gz"
                    )
                    tissue["cov"] = os.path.join(
                        gtex_cov_path, f"{tissue_name}.v8.covariates.txt"
                    )
                else:
                    tissue["exp"] = os.path.join(folder_path, tissue["file"])
            A, gene_count = create_corr_matrix(tissues, session_id)
            n = len(A) / len(tissues)
            if measure == "query":
                query_tissue = request.form.get("query-tissue")
                query_file = request.files.get("query-file")
                query_index = int(request.form.get("tissue-index"))
                s = int(query_index * n)
                e = int(s + n)
                query_gene_list = pd.read_csv(query_file, header=None)[0].tolist()
                query_gene_list = [
                    f"{gene}.{query_index+1}" for gene in query_gene_list
                ]
                genes = A.columns[s:e]
                common_target_genes = np.intersect1d(genes, query_gene_list)
                genes_indices = [
                    i for i, e in enumerate(A.columns) if e in common_target_genes
                ]
                _, g = right_target_global_centrality_t(
                    A.values,
                    num_layers=len(tissues),
                    target_tissue=query_index,
                    target_gene_indices=genes_indices,
                    start=s,
                    end=e,
                    p=0.9,
                )
                result["Centrality"] = g
            elif measure == "global":
                result["Centrality"] = global_centrality(A.values, len(tissues), p=0.9)
            else:
                result["Centrality"] = local_centrality(A.values, len(tissues), p=0.9)
            name_list = [
                [tissues[i]["name"] for _ in range(gene_count[i])]
                for i in range(len(tissues))
            ]
            result["Tissue Name"] = list(itertools.chain(*name_list))
            result["Gene Name"] = [re.sub(r"\.\d+$", "", s) for s in A.columns]
            if measure == "local":
                ranking = pd.DataFrame()
                for tissue in tissues:
                    df = result.loc[result["Tissue Name"] == tissue["name"]]
                    df = df.sort_values(by="Centrality", ascending=False)
                    df["Rank"] = [i + 1 for i in range(len(df))]
                    ranking = pd.concat([ranking, df], ignore_index=True)
                result = ranking
            else:
                result = result.sort_values(by="Centrality", ascending=False)
                result["Rank"] = [i + 1 for i in range(len(result))]
            result = result[["Tissue Name", "Gene Name", "Centrality", "Rank"]]
        except Exception as e:
            print(str(e))

        # print(result)
        result_csv = os.path.join("static", "downloads", session_id + ".csv")
        result.to_csv(result_csv)
        result["Centrality"] = np.round(result["Centrality"], 4)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time / 60, 2)
        # unit = 'sec(s)' if execution_time<1 else 'min(s)'
        # time_taken = f'{execution_time} {unit}'
        time_taken = f"{execution_time} min(s)"
        print("Execution time:", time_taken, "\n")

        return render_template(
            "ranking.html",
            columns=result.columns,
            rows=list(result.values.tolist()),
            path=result_csv,
            time_taken=time_taken,
        )


if __name__ == "__main__":
    #    app.run(debug=True)
    app.run(host="0.0.0.0", port=5000, debug=True)
