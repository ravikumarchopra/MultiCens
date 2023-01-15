from flask import Flask, render_template, request, url_for, redirect, flash
from werkzeug.utils import secure_filename
from util.ranking import getRanking
from util.util import allowed_file
from pathlib import Path
import pandas as pd
import uuid
import os

tissue_list = ['adipose', 'liver', 'muscle', 'intestine', 'thyroid', 'ovary', 'breast', 'hypothalamus', 'kidney', 'pancreas', 'pituitary', 'adrenal', 'uterus', 'stomach', 'heart']
hormone_list = ['adrenaline', 'aldosterone', 'angiotensin', 'cortisol', 'estradiol', 'glucagon', 'insulin', 'leptin', 'norepinephrine', 'progesterone', 'somatotrophin', 'thyroxin', 'vitamind']

tissue_list.sort()
hormone_list.sort()

app = Flask(__name__)
@app.route('/')
def index():
    return render_template('home.html', tissue_list=tissue_list, hormone_list=hormone_list)

@app.route('/get_ranking', methods=['POST'])
def get_ranking():
    if request.method == 'POST':
        if len(request.files)<=0:
            flash('No selected file')
            return redirect(request.url)

        file1 = request.files['input1']
        file2 = request.files['input2']
        source_tissue = request.form.get('source_tissue')
        target_tissue = request.form.get('target_tissue')
        hormone = request.form.get('hormone')
        
        source_path = 'static/uploads/' + uuid.uuid4().hex + Path(file1.filename).suffix
        target_path = 'static/uploads/' + uuid.uuid4().hex + Path(file2.filename).suffix
        file1.save(source_path)
        file2.save(target_path)

        # result = pd.read_csv(source)
        result = pd.DataFrame()
        
        try:
            result = getRanking(source_path, target_path, hormone, source_tissue, target_tissue)
        except:
            pass

        # delete file from the server
        os.remove(source_path)
        os.remove(target_path)
        
        result_csv = 'static/downloads/' + uuid.uuid4().hex + '.csv'
        
        result.to_csv(result_csv, index=False)

        return render_template('ranking.html', columns=result.columns, rows=list(result.values.tolist()), path=result_csv)
    
@app.route('/run_sample', methods=['GET'])
def run_sample():
    if request.method == 'GET':
        source_tissue = 'pancreas'
        target_tissue = 'muscle'
        hormone = 'insulin'
        
        source_path = 'static/downloads/insulin_source_genes_full.csv'
        target_path = 'static/downloads/insulin_target_genes_full.csv'
        
        result = pd.DataFrame()
        
        try:
            result = getRanking(source_path, target_path, hormone, source_tissue, target_tissue)
        except:
            pass
        
        result_csv = 'static/downloads/' + uuid.uuid4().hex + '.csv'
        
        result.to_csv(result_csv, index=False)

        return render_template('ranking.html', columns=result.columns, rows=list(result.values.tolist()), path=result_csv)

if __name__=='__main__':
#    app.run(debug=True)
   app.run(host='0.0.0.0', port=5000, debug=False)

   