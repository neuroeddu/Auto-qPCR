import pandas
import numpy as np
import re


def process(data, colnames, target_sorter, sample_sorter):
	outlier_data = data[data['Outliers'].eq(True)]
	data = data[data['Outliers'].eq(False)]
	# Calculate Mean (Endogenous Control Mean) and SSD for all Controls
	data['NormQuant'] = 0

	control_filter = (data['Control'].eq(True))
	# print("Endogenous Control Quantity Means and SSD")
	# print(data_controls_grouped)

	data_controls_quantity = data[control_filter].groupby(['Sample Name'], sort=False).agg({'Quantity': 'mean'})
	# print("Combined Endogenous Control Quantity Means and SSD")
	# print(data_controls_quantity)

	# Create Normalized Quantity column
	for i , row in enumerate(data_controls_quantity.itertuples(name=None) , 1):
		name_filter = (data['Sample Name'] == row[0]) #control filter used to be here
		for j in data[name_filter].index:
			data.loc[j , 'NormQuant'] = data.loc[j , 'Quantity'] / row[1]

	# Calculate the SEM for technical replicate groups
	targets = data['Target Name'].drop_duplicates(keep='first').values
	mean_sem_result = {}
	for target in targets:
		mean_sem_result[target] = {}
		samples = data[data['Target Name'] == target]['Sample Name'].drop_duplicates(keep='first').values
		for sample in samples:
			target_sample_data = data.loc[(data['Target Name'] == target) & (data['Sample Name'] == sample)]
			mean = target_sample_data['NormQuant'].mean()
			sdt_dev = target_sample_data['NormQuant'].std()
			std_err = target_sample_data['NormQuant'].sem()
			mean_sem_result[target][sample] = (mean , sdt_dev , std_err)
	for i_row , row in data.iterrows():
		if data.at[i_row , 'Sample Name'] in samples and data.at[i_row , 'Sample Name'] in samples and data.at[
			i_row , 'Target Name'] in mean_sem_result and data.at[i_row , 'Sample Name'] in mean_sem_result[
			data.at[i_row , 'Target Name']]:
			data.at[i_row , 'NormMean'] = \
				mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][0]
			data.at[i_row , 'NormSD'] = mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][
				1]
			data.at[i_row , 'NormSEM'] = \
				mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][2]

	# Making the intermediate dataframe
	data = data.append(outlier_data) # add outliers
	# Sorting data
	data = data_sorter(data, target_sorter, sample_sorter)

	cnames = [c.strip().lower() for c in colnames.split(',')]
	clist = []
	for c in data.columns.values.tolist():
		if c.lower() in cnames:
			clist.append(c)
	cols = ['Target Name', 'Sample Name', 'NormQuant', 'NormMean', 'NormSD', 'NormSEM', 'Outliers'] + clist
	df = pandas.DataFrame(columns=cols)
	for item in cols:
		df[item] = data[item]
	df.reset_index(drop=True, inplace=True)

	# remove outliers in summary data
	data = data[data['Outliers'].eq(False)]
	data_output_summary = data.groupby(['Target Name', 'Sample Name'], sort=False).agg(
		{'NormQuant': [np.size , 'mean' , 'std'] , 'NormSEM': 'mean'})

	data_output_summary_w_group = data.groupby(['Target Name', 'Sample Name']+clist, sort=False).agg(
		{'NormQuant': [np.size , 'mean' , 'std'] , 'NormSEM': 'mean'})

	targets = data['Target Name'].drop_duplicates(keep='first').values
	samples = data['Sample Name'].drop_duplicates(keep='first').values

	return df, data_output_summary, data_output_summary_w_group, targets, samples


def data_sorter(data, target_sorter, sample_sorter):
	# define sorter for target name order based on list
	targets = data['Target Name'].drop_duplicates(keep='first').values
	if target_sorter != '':
		targets = [sorter.strip() for sorter in target_sorter.split(',')]
	# remove nan from list
	targets = [t for t in targets if type(t) is not float]
	# Add sorter to dataframe to order Target Names in output files
	sorter_index = dict(zip([g.lower() for g in targets], range(len(targets))))
	data['Target Order'] = data['Target Name'].str.lower().map(sorter_index)
	data.sort_values(['Target Order'], inplace=True)

	# define sorter for sample name order based on list
	sorter = data['Sample Name'].drop_duplicates(keep='first').values
	if sample_sorter != '':
		sorter = [sorter.strip() for sorter in sample_sorter.split(',')]
	# remove nan from list
	sorter = [s for s in sorter if type(s) is not float]
	# Add sorter to dataframe to order Sample Names in output files
	sorter_index = dict(zip([s.lower() for s in sorter], range(len(sorter))))
	data['Sample Name Key'] = data['Sample Name'].str.extract(
			re.compile('(' + '|'.join(sorter) + ')', re.IGNORECASE), expand=False).fillna('')
	data['Sample Order'] = data['Sample Name Key'].str.lower().map(sorter_index)
	data.sort_values(by=['Target Order', 'Sample Order'], inplace=True)

	return data
