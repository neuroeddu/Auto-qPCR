# ID Lab qPCR Analysis
# Zero Clause BSD Copyright (c) 2019 by Iveta Demirova

VERSION = "0.1.7"
QUALITY = ""

import pandas
import pingouin as pg
import numpy as np
from pingouin import pairwise_ttests, multicomp, ttest
import absolute, relative, stability

def process_data(data , model , cgenes , cutoff , max_outliers , csample=None):
    """This filters the data and processes the selected model, returning a list of output dataframes"""

    # Transforms certain columns from string to numeric
    cols = ['CT' , 'Quantity']
    data[cols] = data[cols].apply(pandas.to_numeric , errors='coerce')

    # Marks the Control Genes in a new column in the dataframe
    data['Control'] = data['Target Name'].apply(lambda x: True if str(x) in cgenes else False)

    # Create column 'Ignore' in dataframe to mark rows with NaN values in certain columns
    data['Ignore'] = False
    data['Outliers'] = False
    cols = ['Sample Name' , 'Target Name' , 'Task' , 'Reporter' , 'CT']
    for col in cols:
        data.loc[data[col].isnull() , 'Ignore'] = True

    # Calls the different processing models depending on the model argument
    if model == 'absolute':
        data = cleanup_outliers(data , "Quantity" , cutoff , max_outliers)
        data, data_summary, targets, samples = absolute.process(data)

    elif model == 'relative':
        data = cleanup_outliers(data , "CT" , cutoff , max_outliers)
        data, data_summary, targets, samples  = relative.process(data)

    elif model == 'stability':
        data = cleanup_outliers(data , "CT" , cutoff , max_outliers)
        data, data_summary, targets, samples = stability.process(data , csample)

    return data, data_summary, targets, samples


def cleanup_outliers(d , feature , cutoff , max_outliers):
    """Function to remove outliers based on cutoff and maximum number of outliers,
    by removing the furthest data point in each group when the standard deviation
    is higher than the cutoff"""

    # Calculate SSD for all sample groups
    f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN')
    d1 = d[f].groupby(['Sample Name' , 'Target Name']).agg({'CT': ['std']})
    f = (d1['CT']['std'] > cutoff)
    d2 = d1[f]
    if not d2.empty:
        # Mark all outliers
        for i , row in enumerate(d2.itertuples(name=None) , 1):
            f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
                & (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
            dx_idx = d[f].index
            group_size = len(dx_idx)
            min_size = round(group_size * (1 - max_outliers))
            size = group_size
            if min_size < 2:
                min_size = 2
                h.log(5 , 'Minimum group size must be equal or greater than 2')
            while True:
                f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
                    & (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
                dx = d[f].copy()
                dxg = d[f].groupby(['Sample Name' , 'Target Name']).agg({feature: [np.size , 'std' , 'mean']})
                if dxg[feature]['std'].iloc[0] <= cutoff:
                    # CT std is under the threshold
                    break
                # Will ignore one or all measurements
                size -= 1
                if size < min_size:
                    # Ignore the entire group of measurements
                    # for j in dx_idx:
                    #    d['Ignore'].loc[j] = True
                    break
                # Will remove the measurement which is furthest from the mean
                dx['Distance'] = (dx[feature] - dxg[feature]['mean'].iloc[0]) ** 2
                j = dx.sort_values(by='Distance' , ascending=False).index[0]
                d['Outliers'].loc[j] = True
                d['Ignore'].loc[j] = True

    return (d[(d['Ignore'].eq(False))])

