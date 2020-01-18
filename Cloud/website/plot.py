import plotly.graph_objs as go
from plotly.graph_objs import *
import pandas as pd
import AUTOqPCR
import plotly.offline as py

def main():
    # data = pd.read_csv('/Users/admin/Documents/GitHub/Auto-qPCR/Manuscript/DataSets/stability/STABILITY_TEST.csv',
    #                    skip_blank_lines=True ,
    #                    skipinitialspace=True ,
    #                    engine='python' ,
    #                    encoding="utf-8" ,
    #                    header=46)
    # data1, summary_data, targets, samples = AUTOqPCR.process_data(data,'stability', 'CHR4',  0.3, 0.5, csample='A')
    # fig = plot_by_targets(summary_data, 'stability', targets, samples)

    data = pd.read_csv('/Users/admin/Documents/GitHub/Auto-qPCR/Manuscript/DataSets/Absolute/Gilles_data_combined.csv',
                                          skip_blank_lines=True ,
                                          skipinitialspace=True ,
                                          engine='python' ,
                                          encoding="utf-8" ,
                                          header=46)
    data1, summary_data, targets, samples = AUTOqPCR.process_data(data, 'absolute', 'GAPDH ACTB', 0.3, 0.5)
    # data = pd.read_csv('intermediate_absolute.csv')
    # fig = plot_by_groups(data, 'absolute')

    # data = pd.read_csv('/Users/admin/Downloads/output_absolute.csv')
    fig = plot_by_targets(summary_data, 'absolute', targets, samples)
    fig.show()


def plot_by_targets(dataframe, model, targets, samples):
    data_list = []
    # plot x-axis in desired order
    layout = go.Layout(xaxis=dict(type='category'))

    if model == 'absolute':
        for item in targets:
            data = go.Bar(
                name=item,
                x=list(samples),
                y=dataframe.loc[item, 'NormQuant']['mean'],
                error_y=dict(type='data', array=dataframe.loc[item, 'NormSEM']['mean'])
            )
            data_list.append(data)

    else:
        for item in targets:
            data = go.Bar(
                name=item,
                x=list(samples),
                y=dataframe.loc[item, 'rq']['mean'] ,
                error_y=dict(type='data', array=dataframe.loc[item, 'rqSEM']['mean'])
            )
            data_list.append(data)

        for item in samples:
            data = go.Bar(
                name=item ,
                x=list(targets) ,
                y=dataframe.loc[(slice(None) , item) , 'rq']['mean'] ,
                error_y=dict(type='data' , array=dataframe.loc[(slice(None) , item) , 'rqSEM']['mean'])
            )
            data_list.append(data)

    return data_list, layout


def plot_by_groups(df, model, groups=None, targets=None):
    if groups is None:
        groups = df['Group'].drop_duplicates(keep='first').values.tolist()
    else:
        groups=groups
    if targets is None:
        targets = df['Target Name'].drop_duplicates(keep='first').values.tolist()

    data_list=[]
    layout = go.Layout(xaxis=dict(type='category'))
    if model == 'absolute':
        count = 0
        for t in targets:
            y = []
            st_err = []
            count += 1
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['NormMean'].mean())
                st_err.append(sample['NormMean'].sem())
        data = go.Bar(name=t , x=groups , y=y , error_y=dict(type='data' , array=st_err))
        data_list.append(data)

    else:
        count = 0
        for t in targets:
            y = []
            st_err = []
            count += 1
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['rqMean'].mean())
                st_err.append(sample['rqMean'].sem())

            data=go.Bar(name=t, x=groups, y=y, error_y=dict(type='data', array=st_err))
            data_list.append(data)

    return data_list, layout


if __name__=='__main__':
    main()