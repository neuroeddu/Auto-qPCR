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
	
	control_names = cfg['MODEL']['ControlGenes'].split(', ')

	target_names = list(dict.fromkeys(d['Target Name']))
	targets = dict((target,[]) for target in target_names)

	for control in control_names:
		targets.pop(control, None)

	controls = dict((control,[]) for control in control_names)

	for i, well in enumerate(d['Well']):

		if d['CT'][i] != "Undetermined":

			t = d['Target Name'][i]

			if t in control_names:
				controls[t].append(float(d['CT'][i]))
			else:
				targets[t].append(float(d['CT'][i]))

	control_means = []

	for control, array in controls.items():
		if len(array) > 0:
			control_means.append(np.mean(array))

	control_final_mean = np.mean(control_means)

	for target, array in targets.items():
		for i, item in enumerate(array):
			ct = item 
			delta_ct = ct - control_final_mean
			RQ = np.power(2, -delta_ct)
			targets[target][i] = RQ

	target_means = {}
	target_stderr = {}

	for target, array in targets.items():
		target_means[target] = np.mean(array)
		target_stderr[target] = np.std(array) / np.sqrt(len(array))


	print(controls)

	print(control_means)

	print(control_final_mean)

	print(targets)

	print(target_means)


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