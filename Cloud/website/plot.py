import pandas as pd
import AUTOqPCR
import matplotlib.pyplot as plt
import numpy as np


def plots(dataframe, model, targets, samples, cgenes):
    plots=[]
    # Absolute: bar plots for each gene, a grouped bar plot by genes and a grouped bar plot by cell lines
    if model == 'absolute':
        plot_by_genes = plt.figure(figsize=(40, 20))
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
            # remove endogenous control from individual plots
            if item.lower() not in cgenes.lower().split(','):
                plots.append(plot)
            # grouped plot
            plt.bar(x2 , sample , yerr=list(dataframe.loc[item , 'NormSEM']['mean']) , align='center' ,
                    error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , width=barwidth , edgecolor = 'white', label=item)
            counter += 1
        plt.xticks([i*len(targets) + barwidth*counter/2 for i in range(len(samples))] , samples , rotation='vertical' , fontsize='20')
        plt.xlabel('Targets', fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close(plot_by_genes)
        plots.append(plot_by_genes)

        plot_by_samples = plt.figure(figsize=(40, 20))
        barwidth = 0.75
        counter = 0
        for item in samples:
            target = list(dataframe.loc[(slice(None) , item) , 'NormQuant']['mean'])
            x = np.arange(len(target)) * len(samples) + barwidth * counter
            plt.bar(x , target , yerr=list(dataframe.loc[(slice(None) , item) , 'NormSEM']['mean']) , align='center' ,
                    error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , width=barwidth , edgecolor='white' , label=item)
            counter += 1

        plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))] , targets ,
                   rotation='vertical' , fontsize='20')
        plt.xlabel('Samples', fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_samples)

    # Genomic stability: grouped by chromosomes and by cell lines
    elif model == 'stability2':
        # plot grouped by chromosomes
        plot_by_genes = plt.figure(figsize=(40, 20))
        barwidth = 0.75
        counter = 0
        for item in targets:
            sample = list(dataframe.loc[item , 'rq']['mean'])
            x = np.arange(len(sample)) * len(targets) + barwidth * counter
            plt.bar(x , sample , yerr=list(dataframe.loc[item , 'rqSEM']['mean']) , align='center' ,
                    error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , width=barwidth , edgecolor='white' , label=item)
            counter += 1

        plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))] , samples ,
                   rotation='vertical' , fontsize='20')
        plt.xlabel('Targets' , fontsize='20' , fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_genes)
        # plot grouped by samples
        plot_by_samples = plt.figure(figsize=(40 , 20))
        barwidth = 0.75
        counter = 0
        for item in samples:
            target = list(dataframe.loc[(slice(None) , item) , 'rq']['mean'])
            x = np.arange(len(target)) * len(samples) + barwidth * counter
            plt.bar(x , target , yerr=list(dataframe.loc[(slice(None) , item) , 'rqSEM']['mean']) , align='center' ,
                    error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , width=barwidth , edgecolor='white' , label=item)
            counter += 1

        plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))] , targets ,
                   rotation='vertical' , fontsize='20')
        plt.xlabel('Samples', fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_samples)

    # relative deltaCT and delta delta CT: grouped by genes
    else:
        grouped_plot = plt.figure(figsize=(40 , 20))
        barwidth = 0.75
        counter = 0
        for item in targets:
            plot = plt.figure(figsize=(20 , 20))
            sample = list(dataframe.loc[item , 'rq']['mean'])
            x = np.arange(len(sample))
            x2 = x * len(targets) + barwidth * counter
            plt.bar(x , sample , yerr=list(dataframe.loc[item , 'rqSEM']['mean']) , align='center' ,
                    error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , width=barwidth , edgecolor='white' , label=item)
            plt.xlabel(item , fontweight='bold')
            plt.xticks([i for i in range(len(samples))] , samples , rotation='vertical' , fontsize='20')
            plt.legend(fontsize='20')
            plt.close(plot)
            # remove endogenous control from individual plots
            if item.lower() not in cgenes.lower().split(','):
                plots.append(plot)
            # grouped plot
            plt.bar(x2 , sample , yerr=list(dataframe.loc[item , 'rqSEM']['mean']) , align='center' ,
                    error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , width=barwidth , edgecolor='white' , label=item)
            counter += 1
        plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))], samples,
                   rotation='vertical' , fontsize='20')
        plt.xlabel('Targets', fontsize='20' , fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(grouped_plot)

    return plots


# plot by user defined groups
def plot_by_groups(df, model, targets):
    groups = df['Group'].drop_duplicates(keep='first').values.tolist()

    plots = []
    if model == 'absolute':
        barwidth = 0.75
        # grouped by groups on the x-axis
        counter = 0
        plot_by_group = plt.figure(figsize=(20,20))
        for t in targets:
            y = []
            st_err = []
            x = np.arange(len(groups))*len(targets)+barwidth*counter
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['NormMean'].mean())
                st_err.append(sample['NormMean'].sem())
            plt.bar(x, y, yerr=st_err, error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9), align='center', width=barwidth, edgecolor='white', label=t)

            counter += 1
        plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(groups))] , groups ,
                   rotation='vertical' , fontsize='20')
        plt.xlabel('Groups' , fontsize='20' , fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_group)
        # grouped by targets on the x-axis
        counter = 0
        plot_by_target = plt.figure(figsize=(20 , 20))
        for g in groups:
            y = []
            st_err = []
            x = np.arange(len(targets))*len(groups)+barwidth*counter
            for t in targets:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['NormMean'].mean())
                st_err.append(sample['NormMean'].sem())
            plt.bar(x, y, yerr=st_err, error_kw=dict(lw=0.9, capsize=2, capthick=0.9), align='center', width=barwidth, edgecolor='white', label=g)

            counter += 1
        plt.xticks([i * len(groups) + barwidth * counter / 2 for i in range(len(targets))], targets,
                   rotation='vertical', fontsize='20')
        plt.xlabel('Targets', fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_target)
    else:
        barwidth = 0.75
        # grouped by groups on the x-axis
        counter = 0
        plot_by_group = plt.figure(figsize=(20, 20))
        for t in targets:
            y = []
            st_err = []
            x = np.arange(len(groups))*len(targets)+barwidth*counter
            for g in groups:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['rqMean'].mean())
                st_err.append(sample['rqMean'].sem())
            plt.bar(x , y , yerr=st_err , error_kw=dict(lw=0.9, capsize=2, capthick=0.9), align='center', width=barwidth, edgecolor='white', label=t)

            counter += 1
        plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(groups))], groups,
                   rotation='vertical' , fontsize='20')
        plt.xlabel('Groups' , fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_group)
        # grouped by groups on the x-axis
        counter = 0
        plot_by_target = plt.figure(figsize=(20, 20))
        for g in groups:
            y = []
            st_err = []
            x = np.arange(len(targets)) * len(groups) + barwidth * counter
            for t in targets:
                sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
                y.append(sample['rqMean'].mean())
                st_err.append(sample['rqMean'].sem())
            plt.bar(x , y , yerr=st_err , error_kw=dict(lw=0.9 , capsize=2 , capthick=0.9) , align='center' ,
                    width=barwidth , edgecolor='white' , label=g)

            counter += 1
        plt.xticks([i * len(groups) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='vertical', fontsize='20')
        plt.xlabel('Targets', fontsize='20', fontweight='bold')
        plt.legend(fontsize='20')
        plt.close()
        plots.append(plot_by_target)

    return plots


if __name__=='__main__':
    main()