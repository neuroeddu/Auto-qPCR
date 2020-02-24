from flask import Flask, request, make_response, render_template
import io
import pandas as pd
import AUTOqPCR
import plot
import statistics
from zipfile import ZipFile


app: Flask = Flask(__name__)


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
        data = data.append(filedata , ignore_index=True , sort=True)
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

    data1, summary_data, targets, samples = AUTOqPCR.process_data(data , model , cgenes , cutoff , max_outliers , sample_sorter , csample)
    # , sorter

    # taking lists of samples, targets and groups in the order user want to plot
    otargets = request.form['otargets'].split()
    if len(otargets) != 0:
        targets = otargets

    # making stats csv
    if qty is not None:
        if request.form['option4'] != 'False':
            gcol = data[request.form['gcol']]
            data1['Group'] = gcol
        else:
            groups = request.form['glist'].split()
            data1 = statistics.add_groups(data1, groups)

        anova_dfs , posthoc_dfs = statistics.stats(model, qty , data1, targets , rm , posthoc)
        print(anova_dfs)
        print(posthoc_dfs)
        anova_output = anova_dfs.to_csv(index=False)
        posthoc_output = posthoc_dfs.to_csv(index=False)

    # online plotly plots
    # if sorter is not None:
    #     plots = plot.plot_by_targets(summary_data, model, targets, sorter)
    # else:
    plots = plot.plot_by_targets(summary_data , model , targets , samples)

    # making summary data csv
    output = summary_data.to_csv()
    clean_output = data1.to_csv()

    outfile = io.BytesIO()
    with ZipFile(outfile, 'w') as myzip:
        if qty is not None:
            myzip.writestr('anova_result.csv', anova_output)
            myzip.writestr(posthoc+'_result.csv' , posthoc_output)
        myzip.writestr('clean_data.csv' , clean_output)
        myzip.writestr('summary_data.csv', output)
        for i in range(len(plots)):
            buf = io.BytesIO()
            plots[i].savefig(buf)
            image_name = 'image {}.png'.format(i+1)
            myzip.writestr(image_name, buf.getvalue())
        myzip.close()

    response = make_response(outfile.getvalue())
    response.headers['Content-Type'] = 'application/actet-stream'
    response.headers['Content-Disposition'] = 'attachment; filename=outputs_'+model+'.zip'

    return response


if __name__=='__main__':
    app.debug = True
    app.run(host = '0.0.0.0', port=5000)

