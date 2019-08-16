from flask import Flask, request, make_response, render_template
import io
import pandas as pd
import AUTOqPCR
import json
import plotly
import plotly.graph_objs as go


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
    max_outliers = request.form.get('max_outliers' , type=float)
    csample = request.form['csample']

    result, targets, samples = AUTOqPCR.process_data(data , model , cgenes , cutoff , max_outliers , csample)

    output = result.to_csv()

    response = make_response(output)
    response.headers['Content-Disposition'] = 'attachment; filename=output.csv'

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

        fig.show()
    else:
        fig = go.Figure()
        for item in targets:
            fig.add_trace(go.Bar(
                name=item,
                x=list(samples),
                y=dataframe.loc[item, 'rq']['mean'] ,
                error_y=dict(type='data', array=dataframe.loc[item, 'rqSEM']['mean'])
            ))

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    return graphJSON


if __name__=='__main__':
    app.debug = True
    app.run(host = '0.0.0.0', port=5000)

