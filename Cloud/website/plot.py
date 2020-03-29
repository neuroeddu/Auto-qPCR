import pandas as pd
import AUTOqPCR
import matplotlib.pyplot as plt
import numpy as np


def main():
    data = pd.read_csv('Gilles_data_combined.csv',
                                   skip_blank_lines=True ,
                                   skipinitialspace=True ,
                                   engine='python' ,
                                   encoding="utf-8" ,
                                   header=46)
    data1, summary_data, targets, samples = AUTOqPCR.process_data(data, 'absolute', 'GAPDH ACTB', 0.3, 0.5)
    # targets = ['ACTB', 'GAPDH', 'GRIA1', 'GRIN1', 'KCNJ6', 'SYP']
    # print(df.loc['ACTB', 'NormQuant']['mean'])
    # samples = ['AiW001-2', 'AJC001-5', '522-266-2', 'AJG001C4', 'NCRM1']
    plots = plot_by_targets(summary_data, 'absolute', targets, samples)
    for i in range(len(plots)):
        plots[i].savefig('image {}.jpg'.format(i+1))


def plot_by_targets(dataframe, model, targets, samples):

    plots=[]
    if model == 'absolute':
        for item in targets:
            plot = plt.figure(figsize=(20, 20))
            sample = list(dataframe.loc[item, 'NormQuant']['mean'])
            x = np.arange(len(sample))
            plt.bar(x, sample, yerr=list(dataframe.loc[item, 'NormSEM']['mean']), align='center',
                   error_kw=dict(lw=0.9, capsize=2, capthick=0.9), width=0.75, label=item)

            plt.xlabel(item, fontweight='bold')
            plt.xticks([i for i in range(len(samples))], samples, rotation='vertical', fontsize='20')
            plt.legend(fontsize='20')
            plt.close()
            plots.append(plot)
    else:
        for item in targets:
            plot = plt.figure(figsize=(20, 20))
            sample = list(dataframe.loc[item , 'rq']['mean'])
            x = np.arange(len(sample))
            plt.bar([i for i in x] , sample , yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
                   error_kw=dict(lw=0.9, capsize=2, capthick=0.9), width=0.75, label=item)

            plt.xlabel(item , fontweight='bold')
            plt.xticks([i for i in range(len(samples))] , samples , rotation='vertical' , fontsize='20')
            plt.legend(fontsize='20')

            plt.close()
            plots.append(plot)

        for item in samples:
            plot = plt.figure(figsize=(20, 20))
            target = list(dataframe.loc[(slice(None), item) , 'rq']['mean'])
            x = np.arange(len(target))
            plt.bar([i for i in x] , target , yerr=list(dataframe.loc[(slice(None), item), 'rqSEM']['mean']),
                    align='center', error_kw=dict(lw=0.9, capsize=2, capthick=0.9), width=0.75 , label=item)

            plt.xlabel(item , fontweight='bold')
            plt.xticks([i for i in range(len(targets))] , targets , rotation='vertical' , fontsize='20')
            plt.legend(fontsize='20')

            plt.close()
            plots.append(plot)

    return plots


def plot_by_groups(df, model, groups=None, targets=None):
    if groups is None:
        groups = df['Group'].drop_duplicates(keep='first').values.tolist()
    else:
        groups=groups
    if targets is None:
        targets = df['Target Name'].drop_duplicates(keep='first').values.tolist()

    if model == 'absolute':
        fig , ax = plt.subplots()
        width = 0
        for t in targets:
            y = []
            st_err = []
            x = np.arange(targets)
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['NormMean'].mean())
                st_err.append(sample['NormMean'].sem())

                ax.bar([i + width for i in x] , y , yerr=st_err, align='center')

    else:
        fig , ax = plt.subplots()
        width = 0
        for t in targets:
            y = []
            st_err = []
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['rqMean'].mean())
                st_err.append(sample['rqMean'].sem())

                ax.bar([i + width for i in x] , y , yerr=st_err , align='center')


if __name__=='__main__':
    main()