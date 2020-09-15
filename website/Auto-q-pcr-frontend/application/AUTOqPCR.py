# ID Lab qPCR Analysis
# Zero Clause BSD Copyright (c) 2019 by Iveta Demirova

VERSION = "0.1.7"
QUALITY = ""

import pandas
import numpy as np
from application import absolute , relative , stability
# from tabulate import tabulate


def process_data(data , model , quencher, task, cgenes , cutoff , max_outliers , preservevar, target_sorter=None , sample_sorter=None , csample=None, colnames=None):
	"""This filters the data and processes the selected model, returning a list of output dataframes"""
	# Transforms certain columns from string to numeric
	cols = ['CT' , 'Quantity']
	data[cols] = data[cols].apply(pandas.to_numeric , errors='coerce')

	# Marks the Control Genes in a new column in the dataframe
	data['Control'] = data['Target Name'].apply(lambda x: True if str(x).lower() in cgenes.lower() else False)
	# Create column 'Ignore' in dataframe to mark rows with NaN values in certain columns
	data['Ignore'] = False
	data['Outliers'] = False
	cols = ['Sample Name' , 'Target Name' , 'Task' , 'Reporter' , 'CT']
	for col in cols:
		data.loc[data[col].isnull() , 'Ignore'] = True
	# ignore rows if they correspond to signal from the quencher
	data['Ignore'].mask(data['Reporter'].str.lower() == quencher.lower(), True, inplace=True)

	# Calls the different processing models depending on the model argument
	if model == 'absolute':
		data = cleanup_outliers(data, "Quantity", cutoff, max_outliers, preservevar, task)
		data, data_summary, data_summary_w_group, targets, samples = absolute.process(data, colnames, target_sorter, sample_sorter)

	elif model == 'relative_dCT':
		data = cleanup_outliers(data, "CT", cutoff, max_outliers, preservevar, task)
		data, data_summary, data_summary_w_group, targets, samples = relative.process(data, colnames, target_sorter, sample_sorter)

	else:
		data = cleanup_outliers(data, "CT", cutoff, max_outliers, preservevar, task)
		data, data_summary, data_summary_w_group, targets, samples = stability.process(model, data, csample, colnames, target_sorter, sample_sorter)

	return data, data_summary, data_summary_w_group, targets, samples


def cleanup_outliers(d , feature , cutoff , max_outliers, preservevar, task):
	"""Function to remove outliers based on cutoff and maximum number of outliers,
	by removing the furthest data point in each group when the standard deviation
	is higher than the cutoff"""

	# Calculate SSD for all sample groups
	f = (d['Ignore'].eq(False)) & (d['Task'].str.lower() == task.lower())
	d1 = d[f].groupby(['Sample Name' , 'Target Name']).agg({'CT': ['std']})
	f = (d1['CT']['std'] > cutoff)
	d2 = d1[f]
	if not d2.empty:
		# Mark all outliers
		for i , row in enumerate(d2.itertuples(name=None) , 1):
			f = (d['Ignore'].eq(False)) & (d['Task'].str.lower() == task.lower()) \
				& (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
			dx_idx = d[f].index
			group_size = len(dx_idx)
			min_size = round(group_size * (1 - max_outliers))
			size = group_size
			if min_size < 2:
				min_size = 2
			while True:
				f = (d['Ignore'].eq(False)) & (d['Task'].str.lower() == task.lower()) \
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
				# check if the outlier should be kept if mean has high variation
				if preservevar == 'True':
					if abs((dxg[feature]['mean'].iloc[0]-dx[feature].median())/dx[feature].median()) < 0.1:
						d['Outliers'].loc[j] = False

	return d[(d['Ignore'].eq(False))]
