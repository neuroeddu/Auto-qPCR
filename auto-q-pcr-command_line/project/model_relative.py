import os
import re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import helpers as h
import config #as c
import warnings

def run_model(wdir, d, cfg):
    # Clean up column names


    d.columns = [x.strip() for x in d.columns]

    # Add filter columns
    d['Ignore'] = False
    d['Control'] = False
    d['rq'] = np.nan
    d['rqSEM'] = np.nan
    d['rqSD'] = np.nan
    d['rqMeanBio'] = np.nan
    d['rqSEMBio'] = np.nan
    d['rqSDBio'] = np.nan
    
    # define sorter for sample name order based on list in local config file
    # this list is case sensitive
    sorter = [x.strip() for x in cfg['MODEL']['SampleOrder'].split(',')]
    # print(sorter)
    
    # Convert all values in selected columns to numeric
    cols = ['CT']
    d[cols] = d[cols].apply(pd.to_numeric, errors='coerce')
    
    # Mark the data for the controls defined in the config file in a new column
    # This will allow filtering of the dataframe and different treatment of controls
    ctrls = set([x.strip() for x in cfg['MODEL']['ControlGenes'].split(',')])    
    d['Control'] = d['Target Name'].apply(lambda x: True if x in ctrls else False)
    
    # Create column 'Ignore' in dataframe to mark rows with NaN values in certain columns 
    cols = ['Sample Name', 'Target Name', 'Task', 'Reporter', 'CT']
    for col in cols:
        d.loc[d[col].isnull(), 'Ignore'] = True
    
    # Add a Sample Name Key column for biological sample grouping
    # Combines Sample Names with same name except for numbers at end
    rx = re.compile(r'[ -][0-9]+|[0-9]$', re.IGNORECASE)
    d['Sample Name Key'] = d['Sample Name'].apply(lambda x: rx.sub('', x) if \
    pd.notnull(x) else np.nan)
    
    # Calling function defined at the end of document
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cleanup_ouliers(d, cfg)

    # Add sorter to dataframe to order Sample Names in output files
    targets = set(d['Target Name'])
    for target in targets:
        sorter_index = dict(zip(sorter, range(len(sorter))))
        d['Sample Order'] = d['Sample Name'].map(sorter_index)
        d.sort_values(['Sample Order'], inplace = True)
    
    # Calculate CT Mean (Endogenous Control Mean) and SSD for all Controls
    f1 = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') & (d['Control'].eq(True))
    d1 = d[f1].groupby(['Target Name','Sample Name']).agg({'CT': [np.size, 'mean', 'std']})
    d1_1 = d[f1].groupby(['Target Name','Sample Name']).agg({'CT': [np.size, 'mean', 'std']})
    s = "Endogenous Control CT Means and SSD"
    if h.verbosity == h.LOG_DEBUG:    
        h.bprint(s, 74)
        print(d1_1)
    d2 = d1.groupby(['Sample Name']).agg({('CT', 'mean'):'mean'})
    d2_1 = d1_1.groupby(['Sample Name']).agg({('CT', 'mean'):'mean'})    
    s = "Combined Endogenous Control CT Means and SSD"
    if h.verbosity == h.LOG_DEBUG:
        h.bprint(s, 74)
        print(d2_1)
    
    #Create rq column        
    for i, row in enumerate(d2.itertuples(name = None), 1):
        f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
            & (d['Control'].eq(False)) & (d['Sample Name'] == row[0])
        for j in d[f].index:
            d.loc[j, 'rq'] = np.power(2, -(d.loc[j, 'CT'] - row[1]))
            
    # Calculate the SEM for technical and biological replicate groups
    targets = set(d['Target Name'])
    mean_sem_result = {}
    targets = set(d['Target Name'])
    for target in targets:
        mean_sem_result[target] = {}
        samples = set(d[d['Target Name'] == target]['Sample Name'])
        for sample in samples:
            target_sample_data = d[(d['Target Name'] == target) & (d['Sample Name'] == sample) & d['Ignore'].eq(False)]
            mean = target_sample_data['rq'].mean()
            sdt_dev = target_sample_data['rq'].std()
            std_err = target_sample_data['rq'].sem()
            mean_sem_result[target][sample] = (mean, sdt_dev, std_err)
    for i_row, row in d.iterrows():
        if d.at[i_row, 'Sample Name'] in samples and d.at[i_row, 'Sample Name'] in samples and d.at[i_row, 'Target Name'] in mean_sem_result and d.at[i_row, 'Sample Name'] in mean_sem_result[d.at[i_row, 'Target Name']]:
            d.at[i_row, 'rq'] = mean_sem_result[d.at[i_row, 'Target Name']][d.at[i_row, 'Sample Name']][0]
            d.at[i_row, 'rqSD'] = mean_sem_result[d.at[i_row, 'Target Name']][d.at[i_row, 'Sample Name']][1]
            d.at[i_row, 'rqSEM'] = mean_sem_result[d.at[i_row, 'Target Name']][d.at[i_row, 'Sample Name']][2]
    
    mean_sem_result_bio = {}
    targets = set(d['Target Name'])
    for target in targets:
        mean_sem_result_bio[target] = {}
        samples = set(d[d['Target Name'] == target]['Sample Name Key'])
        for sample in samples:
            target_sample_data = d[(d['Target Name'] == target) & (d['Sample Name Key'] == sample) & d['Ignore'].eq(False)]
            mean = target_sample_data['rq'].mean()
            sdt_dev = target_sample_data['rq'].std()
            std_err = target_sample_data['rq'].sem()
            mean_sem_result_bio[target][sample] = (mean, sdt_dev, std_err)
    for i_row, row in d.iterrows():
        if d.at[i_row, 'Sample Name Key'] in samples and d.at[i_row, 'Sample Name'] in samples and d.at[i_row, 'Target Name'] in mean_sem_result_bio and d.at[i_row, 'Sample Name Key'] in mean_sem_result_bio[d.at[i_row, 'Target Name']]:
            d.at[i_row, 'rqMeanBio'] = mean_sem_result_bio[d.at[i_row, 'Target Name']][d.at[i_row, 'Sample Name Key']][0]
            d.at[i_row, 'rqSDBio'] = mean_sem_result_bio[d.at[i_row, 'Target Name']][d.at[i_row, 'Sample Name Key']][1]
            d.at[i_row, 'rqSEMBio'] = mean_sem_result_bio[d.at[i_row, 'Target Name']][d.at[i_row, 'Sample Name Key']][2]
    
    f2 = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') & (d['Control'].eq(False))
    d3 = d[f2].groupby(['Target Name','Sample Name']).agg({'rq': [np.size, 'mean'], 'rqSD': 'mean', 'rqSEM': 'mean'})
    s = "Mean and SSD for all sample groups"
    
    if h.verbosity == h.LOG_DEBUG:
        h.bprint(s, 85)
        print(d3.to_string())
    fn = Path(wdir).joinpath("technical_group_rq.xlsx")
    d3.to_excel(fn, encoding = cfg['FILE']['Encoding'])
    
    # Calculate Mean and SSD for all sample groups
    f3 = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') & (d['Control'].eq(False))
    d4 = d[f3].groupby(['Target Name','Sample Name Key']).agg({'CT': ['mean'], 'rq': [np.size, 'mean', 'std'], 'rqSEMBio': 'mean'})

    #print(d4.to_string())
    
    #fn = Path(wdir).joinpath("biological_group_rq.xlsx")
    #d4.to_excel(fn, encoding = cfg['FILE']['Encoding'])
    
    # Plots the Normalized Quantities for all targets in one figure
    f5 = (d['Task'] == 'UNKNOWN') & (d['Control'].eq(False))
    
    

    for target_name in set(d[f5]['Target Name']):
        names = []
        values = []
        sem = []
        samples = d[f5][d[f5]['Target Name'] == target_name]
        for sample_name in sorter:
            if sample_name in samples['Sample Name'].values:
                names.append(sample_name)
                try:
                    values.append(samples[samples['Sample Name'] == sample_name]['rq'].iat[0])
                    sem.append(samples[samples['Sample Name'] == sample_name]['rqSEM'].iat[0])
                except Exception:
                    values.append(0)
                    sem.append(0)
        for sample_name in samples['Sample Name'].values:
            if sample_name not in names:
                names.append(sample_name)
                try:
                    values.append(samples[samples['Sample Name'] == sample_name]['rq'].iat[0])
                    sem.append(samples[samples['Sample Name'] == sample_name]['rqSEM'].iat[0])
                except Exception:
                    values.append(0)
                    sem.append(0)
        fig, ax = plt.subplots()
        ax.bar(np.arange(len(values)), values, tick_label=names, yerr=sem, capsize=2)
        ax.set_ylabel('Normalized Expression (+/- SEM)')
        ax.set_xticklabels(names, rotation=90)
        ax.set_title(target_name)
        fig.tight_layout()
        fn = Path(wdir).joinpath(f'{target_name}.png'.replace('/', '-'))
        try: os.remove(fn)
        except Exception: pass
        fig.savefig(fn, bbox_inches='tight')
        if h.verbosity == h.LOG_DEBUG:
            plt.show()
        plt.close(fig)
    
    '''fig = sns.FacetGrid(d[f5], row = 'Target Name', sharex=False, sharey=False, height = 5, aspect = 2)
    fig.map(sns.barplot, 'Sample Name', 'NormQuant', ci=None, order= sorter, errwidth=1, capsize=0.5);
    fig.set_xticklabels(rotation=90)
    fig.set_ylabels('Normalized Expression')
    fig.fig.subplots_adjust(wspace=0, hspace=.95)
    plt.show()
    fig.savefig("graphs.svg")'''

    

#Clean outliers from sample and control groups
def cleanup_ouliers(d, cfg):
    # Calculate SSD for all sample groups
    f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN')
    d1 = d[f].groupby(['Sample Name', 'Target Name']).agg({'CT': ['std']})
    f = (d1['CT']['std'] > cfg['MODEL']['OutlierCutoff'])
    d2 = d1[f]
    if not d2.empty:
        # Mark all outliers
        for i, row in enumerate(d2.itertuples(name = None), 1):
            f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
                & (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
            dx_idx = d[f].index
            group_size = len(dx_idx)
            min_size = round(group_size * (1 - cfg['MODEL']['MaxOutliers']))
            size = group_size
            if min_size < 2:
                min_size = 2
                h.log(5, 'Minimum group size must be equal or greater than 2')
            while True:
                f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
                    & (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
                dx = d[f].copy()
                dxg = d[f].groupby(['Sample Name', 'Target Name']).agg({'CT': [np.size, 'std', 'mean']})
                if dxg['CT']['std'].iloc[0] <= cfg['MODEL']['OutlierCutoff']:
                    # CT std is under the threshold
                    break
                # Will ignore one or all measurements
                size -= 1
                if size < min_size:
                    # Ignore the entire group of measurements
                    #for j in dx_idx:
                    #    d['Ignore'].loc[j] = True
                    break
                # Will remove the measurement which is furthest from the mean
                dx['Distance'] = (dx['CT'] - dxg['CT']['mean'].iloc[0])**2
                j = dx.sort_values(by = 'Distance', ascending = False).index[0]
                d['Ignore'].loc[j] = True
