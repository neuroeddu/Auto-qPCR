from flask import Flask, request, make_response, render_template
import io
import pandas as pd
import AUTOqPCR
from datetime import datetime
import tkinter as tk
from tkinter import ttk
app = Flask(__name__)


@app.route('/')
def form():
    return render_template('layout.html')


@app.route('/download', methods=["POST"])
def transform_view():

    files = request.files.getlist('file[]')
    if len(files) == 0:
        return "No file"
    # Creates empty data frame
    data = pd.DataFrame()

    for item in files:
        stream = io.StringIO(item.stream.read().decode("utf-8"), newline="")
        i = 0
        while True:
            h = stream.readline()
            if not h:
                i = -1
                break
            #if not h.replace(",",""):
                #i += 1
            if h.upper().startswith("WELL"):
                break
            i += 1
        if i == -1:
            print("No header found in '{}'".format(item))
            continue
        else:
            print("Header found at {} in '{}'".format(i, item))
        #print(i)

        stream.seek(0)
        filedata = pd.read_csv(stream ,
                                   skip_blank_lines=True ,
                                   skipinitialspace=True ,
                                   engine='python' ,
                                   encoding="utf-8" ,
                                   header= i)

        # print(filedata)

        data = data.append(filedata, ignore_index=True , sort=True)
        data['filename'] = item.filename
        #stream.seek(0)

    model = request.form['option']
    cgenes = request.form['cgenes']
    cutoff = request.form.get('cutoff' , type=float)
    max_outliers = request.form.get('max_outliers', type=float)
    csample = request.form['csample']
    sample_sorter = request.form['sample_sorter']
    # add plotting option if choose relative
    # plotting_option = request.form['plot']
    qty = request.form.get('quantity', type=int)
    rm = request.form['option2']
    posthoc = request.form['option3']

    data1, summary_data, targets, samples, sorter = AUTOqPCR.process_data(data , model , cgenes , cutoff , max_outliers , sample_sorter , csample)


    # taking lists of samples, targets and groups in the order user wants to plot
    otargets = request.form['otargets'].split()
    if len(otargets) != 0:
        targets = otargets

  # making summary data csv
    output = summary_data.to_csv()
    clean_output = data1.to_csv()
# try to save the files
    
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d-%H-%M-%S")

    response = make_response(clean_output)
    response.headers['Content-Disposition'] = 'attachment; filename= clean_output+' +dt_string+ '.csv'
    return response

if __name__=='__main__':
    app.debug = True
    app.run(host = '0.0.0.0', port=5000)

#
