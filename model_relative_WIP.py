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
	
	controls = cfg['MODEL']['ControlGenes'].split(', ')

	averageControls = 0

	for i, well in enumerate(d['Well']):
		target_name = d['Target Name'][i]
		if d['CT'][i] != 'Undetermined':
			ct = float(d['CT'][i])

			for control in controls:
				if target_name == control:
					averageControls += ct

	averageControls = averageControls / len(d['Well'])

	RQ = []

	for i, well in enumerate(d['Well']):
		target_name = d['Target Name'][i]

		RQ.append('None')

		if d['CT'][i] != 'Undetermined':
			ct = float(d['CT'][i])

			if target_name != control in controls:
				RQ[i] = np.power(2.0, -(ct - averageControls))

	triplicates = list(dict.fromkeys(d['Target Name']))

	targetRQ = dict((triplicate,[]) for triplicate in triplicates)

	for i, well in enumerate(d['Well']):
		target_name = d['Target Name'][i]

		if RQ[i] != 'None':
			if target_name in triplicates:
				targetRQ[target_name].append(int(RQ[i]))


	for key, array in targetRQ.items():
		print(sum(array)/len(array))
		print(np.std(array)/ np.sqrt(3.0))




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
                dxg = d[f].groupby(['Sample Name', 'Target Name']).agg({'CT': [np.size, 'std'], 'Quantity': ['mean']})
                if dxg['CT']['std'].iloc[0] <= cfg['MODEL']['OutlierCutoff']:
                    # CT std is under the threshold
                    break
                # Will ignore one or all measurements
                size -= 1
                if size < min_size:
                    # Ignore the entire group of measurements
                    for j in dx_idx:
                        d['Ignore'].loc[j] = True
                    break
                # Will remove the measurement which is furthest from the mean
                dx['Distance'] = (dx['Quantity'] - dxg['Quantity']['mean'].iloc[0])**2
                j = dx.sort_values(by = 'Distance', ascending = False).index[0]
                d['Ignore'].loc[j] = True