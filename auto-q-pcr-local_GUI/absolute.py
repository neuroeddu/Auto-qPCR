import pandas
import numpy as np


def process(data, colnames=None):

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

	#
	# Making the intermediate dataframe
	data = data.append(outlier_data)
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

	data_output_summary = data.groupby(['Target Name', 'Sample Name'], sort=False).agg(
		{'NormQuant': [np.size , 'mean' , 'std'] , 'NormSEM': 'mean'})

	return df, data_output_summary, targets, samples
