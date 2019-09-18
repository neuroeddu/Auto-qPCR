import plotly.graph_objs as go
from plotly.graph_objs import *
from plotly.subplots import make_subplots
import pandas as pd
import AUTOqPCR


def main():
    # data = pd.read_csv('/Users/admin/Documents/GitHub/Auto-qPCR/Manuscript/DataSets/stability/STABILITY_TEST.csv',
    #                    skip_blank_lines=True ,
    #                    skipinitialspace=True ,
    #                    engine='python' ,
    #                    encoding="utf-8" ,
    #                    header=46)
    # data1, summary_data, targets, samples = AUTOqPCR.process_data(data,'stability', 'CHR4',  0.3, 0.5, csample='A')
    # fig = plot_by_targets(summary_data, 'stability', targets, samples)

    data = pd.read_csv('intermediate_absolute.csv')
    fig = plot_by_groups(data, 'absolute')
    fig.show()


def plot_by_targets(dataframe, model, targets, samples, option=None):
    if model == 'absolute':
        fig = go.Figure()
        for item in targets:
            fig.add_trace(go.Bar(
                name=item,
                x=list(samples),
                y=dataframe.loc[item, 'NormQuant']['mean'],
                error_y=dict(type='data', array=dataframe.loc[item, 'NormSEM']['mean'])
            ))

    elif model == 'relative':
        fig = go.Figure()
        if option == 'targets':
            for item in targets:
                fig.add_trace(go.Bar(
                    name=item,
                    x=list(samples),
                    y=dataframe.loc[item, 'rq']['mean'] ,
                    error_y=dict(type='data', array=dataframe.loc[item, 'rqSEM']['mean'])
                ))
        else:
            for item in samples:
                fig.add_trace(go.Bar(
                    name=item ,
                    x=list(targets) ,
                    y=dataframe.loc[(slice(None) , item) , 'rq']['mean'] ,
                    error_y=dict(type='data' , array=dataframe.loc[(slice(None) , item) , 'rqSEM']['mean'])
                ))
    else:
        fig = go.Figure()
        for item in samples:
            fig.add_trace(go.Bar(
                name=item ,
                x=list(targets) ,
                y=dataframe.loc[(slice(None), item) , 'rq']['mean'] ,
                error_y=dict(type='data' , array=dataframe.loc[(slice(None), item), 'rqSEM']['mean'])
            ))

    fig.update_xaxes(showline=True , linewidth=2 , linecolor='black' , showgrid=False)
    fig.update_yaxes(showline=True , linewidth=2 , linecolor='black' , showgrid=False)

    return fig


def plot_by_groups(df, model):
    groups = df['Group'].drop_duplicates(keep='first').values.tolist()
    targets = df['Target Name'].drop_duplicates(keep='first').values.tolist()
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

            fig = go.Figure(data=go.Bar(name=t, x=groups, y=y, error_y=dict(type='data', array=st_err)))

            # Figure styling
            fig.update_layout(height=600, width=400, title_text=t)
            fig.update_xaxes(showline=True , linewidth=2 , linecolor='black', showgrid=False)
            fig.update_yaxes(showline=True , linewidth=2 , linecolor='black', showgrid=False)

            fig.show()

    # else:
    #     fig = go.Figure()
    #     for item in targets:
    #         fig.add_trace(go.Bar(
    #             name=item ,
    #             x=list(groups) ,
    #             y=df[df['Target Name'].eq(item)]['rq'] ,
    #             error_y=dict(type='data' , array=df[df['Target Name'].eq(item)]['rq'])
    #         ))

    return fig


if __name__=='__main__':
    main()