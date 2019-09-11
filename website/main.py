from flask import Flask, request, make_response, render_template
import io
import pandas as pd
import AUTOqPCR
import plotly.graph_objs as go
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
        stream = io.StringIO(item.stream.read().decode("utf-8") , newline="")
        filedata = pd.read_csv(stream ,
                                   skip_blank_lines=True ,
                                   skipinitialspace=True ,
                                   engine='python' ,
                                   encoding="utf-8" ,
                                   header=46)

        data = data.append(filedata , ignore_index=True , sort=True)
        stream.seek(0)

    model = request.form['option']
    cgenes = request.form['cgenes']
    cutoff = request.form.get('cutoff' , type=float)
    max_outliers = request.form.get('max_outliers', type=float)
    csample = request.form['csample']
    qty = request.form.get('quantity', type=int)
    rm = request.form['option2']
    posthoc = request.form['option3']

    data1, summary_data, targets, samples = AUTOqPCR.process_data(data , model , cgenes , cutoff , max_outliers , csample)
    # making stats csv
    f = open('stats_result.csv', 'a')
    anova_dfs, posthoc_dfs = AUTOqPCR.stats(qty, data1, targets, rm, posthoc)
    for item in anova_dfs:
        item.to_csv(f)
    for item in posthoc_dfs:
        item.to_csv(f)
    f.close()

    fig = create_plot(summary_data, model, targets, samples)
    fig.show()
    # making summary data csv
    output = summary_data.to_csv()
    outfile = io.BytesIO()
    with ZipFile(outfile, 'w') as myzip:
        myzip.writestr('stats_result.csv', f)
        myzip.writestr('summary_data.csv', output)
        myzip.close()

    response = make_response(outfile.getvalue())
    response.headers['Content-Type'] = 'application/actet-stream'
    response.headers['Content-Disposition'] = 'attachment; filename=outputs_'+model+'.zip'

    return response


def create_plot(dataframe, model, targets, samples):
    if model == 'absolute':
        fig = go.Figure()
        for item in targets:
            fig.add_trace(go.Bar(
                name=item,
                x=list(samples),
                y=dataframe.loc[item, 'NormQuant']['mean'],
                error_y=dict(type='data', array=dataframe.loc[item, 'NormSEM']['mean'])
            ))

    else:
        fig = go.Figure()
        for item in targets:
            fig.add_trace(go.Bar(
                name=item,
                x=list(samples),
                y=dataframe.loc[item, 'rq']['mean'] ,
                error_y=dict(type='data', array=dataframe.loc[item, 'rqSEM']['mean'])
            ))

    #graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    return fig


if __name__=='__main__':
    app.debug = True
    app.run(host = '0.0.0.0', port=5000)

