import pandas as pd
import AUTOqPCR
import matplotlib.pyplot as plt
import numpy as np


def plots(dataframe, model, targets, samples):

    plots=[]
    # Absolute: bar plots for each gene and a grouped bar plot by genes
    if model == 'absolute':
        grouped_plot = plt.figure(figsize=(40, 20))
        # barwidth of grouped plot
        barwidth = 0.75
        counter = 0
        for item in targets:
            plot = plt.figure(figsize=(20, 20))
            sample = list(dataframe.loc[item, 'NormQuant']['mean'])
            x = np.arange(len(sample))
            x2 = x*len(targets) + barwidth*counter
            # set color, width, edgecolor, etc.
            plt.bar(x, sample, yerr=list(dataframe.loc[item, 'NormSEM']['mean']), align='center',
                   error_kw=dict(lw=0.9, capsize=2, capthick=0.9), width=barwidth, label=item)
            plt.xlabel(item, fontweight='bold')
            plt.xticks([i for i in range(len(samples))], samples, rotation='vertical', fontsize='20')
            plt.legend(fontsize='20')
            plt.close(plot)
            plots.append(plot)
            # grouped plot
            plt.bar(x2 , sample , yerr=list(dataframe.loc[item , 'NormSEM']['mean']) , align='center' ,
                    error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , width=barwidth , edgecolor = 'white', label=item)
            counter += 1
        plt.xticks([i*len(targets) + barwidth*counter/2 for i in range(len(samples))] , samples , rotation='vertical' , fontsize='20')
        plt.xlabel('Targets',fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close(grouped_plot)
        plots.append(grouped_plot)
    # Genomic stability: grouped by chromosomes
    elif model == 'stability2':
        plot = plt.figure(figsize=(40, 20))
        # barwidth of grouped plot
        barwidth = 0.75
        counter = 0
        for item in targets:
            sample = list(dataframe.loc[item , 'rq']['mean'])
            x = np.arange(len(sample))*len(targets) + barwidth*counter
            plt.bar(x , sample , yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
                   error_kw=dict(lw=0.9, capsize=2, capthick=0.9), width=barwidth, edgecolor='white', label=item)
            counter += 1
        plt.xticks([i*len(targets) + barwidth*counter/2 for i in range(len(samples))], samples, rotation='vertical', fontsize='20')
        plt.xlabel('Targets', fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot)
    # relative deltaCT and delta delta CT: grouped by genes and cell lines
    else:
        # plot grouped by genes
        plot_by_genes = plt.figure(figsize=(40 , 20))
        barwidth = 0.75
        counter = 0
        for item in targets:
            sample = list(dataframe.loc[item , 'rq']['mean'])
            x = np.arange(len(sample))*len(targets) + barwidth*counter
            plt.bar(x , sample , yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
                   error_kw=dict(lw=0.9, capsize=2, capthick=0.9), width=barwidth, edgecolor='white', label=item)
            counter += 1

        plt.xticks([i*len(targets) + barwidth*counter/2 for i in range(len(samples))], samples, rotation='vertical', fontsize='20')
        plt.xlabel('Targets',fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_genes)
        # plot grouped by samples
        plot_by_samples = plt.figure(figsize=(40 , 20))
        barwidth = 0.75
        counter = 0
        for item in samples:
            target = list(dataframe.loc[(slice(None), item) , 'rq']['mean'])
            x = np.arange(len(target))*len(samples) + barwidth*counter
            plt.bar(x , target , yerr=list(dataframe.loc[(slice(None), item), 'rqSEM']['mean']), align='center',
                    error_kw=dict(lw=0.9, capsize=2, capthick=0.9), width=barwidth, edgecolor='white', label=item)
            counter += 1

        plt.xticks([i*len(samples) + barwidth*counter/2 for i in range(len(targets))], targets, rotation='vertical', fontsize='20')
        plt.xlabel('Samples', fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_samples)

    return plots


# If user input groups in the stats section
def plot_by_groups(df, model, groups=None):
    if groups is None:
        groups = df['Group'].drop_duplicates(keep='first').values.tolist()
    else:
        groups=groups
    targets = df['Target Name'].drop_duplicates(keep='first').values.tolist()

    if model == 'absolute':

        fig, ax = plt.subplots()
        width = 0
        for t in targets:
            y = []
            st_err = []
            x = np.arange(len(targets))
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['NormMean'].mean())
                st_err.append(sample['NormMean'].sem())

                ax.bar([i + width for i in x] , y , yerr=st_err, align='center')
    else:
        fig, ax = plt.subplots()
        width = 0
        for t in targets:
            y = []
            st_err = []
            x = np.arange(len(targets))
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['rqMean'].mean())
                st_err.append(sample['rqMean'].sem())

                ax.bar([i + width for i in x] , y , yerr=st_err , align='center')

    return fig


if __name__=='__main__':
    main()