# ID Lab qPCR Analysis
# Zero Clause BSD Copyright (c) 2019 by Iveta Demirova

VERSION = "0.1.5"
QUALITY = ""

from sys import argv
from argparse import ArgumentParser, FileType
import pandas
import numpy as np

def main():

	#Gets arguments
	args = get_cmd_args()

	#Creates empty data frame
	data = pandas.DataFrame()

	#Reads and appends all data from input files to the data dataframe
	for item in args['file']:
		filedata = pandas.read_csv(item,
					skip_blank_lines = True,
					skipinitialspace = True,
					engine = 'python',
					encoding = args['encoding'],
					header = args['header'] - 1)

		data = data.append(filedata, ignore_index = True, sort=True)

	#Runs the main process model with input arguments

	results = process_data(data, args['mod'], args['cgenes'], args['ocutoff'], args['omax'], args['csample'])

	print(results)

	results.to_excel(args['output'] + "output3.xlsx", encoding = args['encoding'])
	#results.to_csv(args['output'] + "output.csv")



def process_data(data, model, cgenes, cutoff, max_outliers, csample = None):
	"""This filters the data and processes the selected model, returning a list of output dataframes"""

	#Transforms certain columns from string to numeric
	cols = ['CT','Quantity']
	data[cols] = data[cols].apply(pandas.to_numeric, errors='coerce')


	#Marks the Control Genes in a new column in the dataframe
	data['Control'] = data['Target Name'].apply(lambda x: True if str(x) in cgenes else False)


	# Create column 'Ignore' in dataframe to mark rows with NaN values in certain columns 
	data['Ignore'] = False
	cols = ['Sample Name', 'Target Name', 'Task', 'Reporter', 'CT']
	for col in cols:
		data.loc[data[col].isnull(), 'Ignore'] = True

	

	# Calls the different processing models depending on the model argument
	if model == 'absolute':
		data = cleanup_outliers(data, "Quantity", cutoff, max_outliers)
		results = process_absolute(data)
	
	elif model == 'relative':
		data = cleanup_outliers(data, "CT", cutoff, max_outliers)
		results = process_relative(data)

	elif model == 'stability':
		data = cleanup_outliers(data, "CT", cutoff, max_outliers)
		results = process_stability(data, csample)

	return results


def process_absolute(data):

	# Calculate Mean (Endogenous Control Mean) and SSD for all Controls

	data['NormQuant'] = 0

	control_filter = (data['Control'].eq(True))
	data_controls_grouped = data[control_filter].groupby(['Target Name','Sample Name']).agg({'Quantity': [np.size, 'mean', 'std']})

	print("Endogenous Control Quantity Means and SSD")
	print(data_controls_grouped)

	data_controls_quantity = data[control_filter].groupby(['Sample Name']).agg({'Quantity':'mean'})

	print("Combined Endogenous Control Quantity Means and SSD")
	print(data_controls_quantity)

	#Create Normalized Quantity column
	for i, row in enumerate(data_controls_quantity.itertuples(name = None), 1):
		name_filter = (data['Control'].eq(False)) & (data['Sample Name'] == row[0])
		for j in data[name_filter].index:
			data.loc[j, 'NormQuant'] = data.loc[j, 'Quantity'] / row[1]

	# Calculate the SEM for technical replicate groups



	targets = set(data['Target Name'])
	mean_sem_result = {}
	for target in targets:
		mean_sem_result[target] = {}
		samples = set(data[data['Target Name'] == target]['Sample Name'])
		for sample in samples:
			target_sample_data = data[(data['Target Name'] == target) & (data['Sample Name'] == sample)]
			mean = target_sample_data['NormQuant'].mean()
			sdt_dev = target_sample_data['NormQuant'].std()
			std_err = target_sample_data['NormQuant'].sem()
			mean_sem_result[target][sample] = (mean, sdt_dev, std_err)
	for i_row, row in data.iterrows():
		if data.at[i_row, 'Sample Name'] in samples and data.at[i_row, 'Sample Name'] in samples and data.at[i_row, 'Target Name'] in mean_sem_result and data.at[i_row, 'Sample Name'] in mean_sem_result[data.at[i_row, 'Target Name']]:
			data.at[i_row, 'NormMean'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][0]
			data.at[i_row, 'NormSD'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][1]
			data.at[i_row, 'NormSEM'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][2]

	data_output_summary = data[(data['Control'].eq(False))].groupby(['Target Name','Sample Name']).agg({'NormQuant': [np.size, 'mean', 'std'], 'NormSEM': 'mean'})

	return data_output_summary




def process_relative(data):
	

	# Calculate CT Mean (Endogenous Control Mean) and SSD for all Controls

	control_filter = (data['Control'].eq(True))
	data_controls_grouped = data[control_filter].groupby(['Target Name','Sample Name']).agg({'CT': [np.size, 'mean', 'std']})

	print("Endogenous Control CT Means and SSD")
	print(data_controls_grouped)

	data_controls_ct = data[control_filter].groupby(['Sample Name']).agg({'CT':'mean'})

	print("Combined Endogenous Control CT Means and SSD")
	print(data_controls_ct)

	#Create rq column
	for i, row in enumerate(data_controls_ct.itertuples(name = None), 1):
		name_filter = (data['Control'].eq(False)) & (data['Sample Name'] == row[0])
		for j in data[name_filter].index:
			data.loc[j, 'rq'] = np.power(2, -(data.loc[j, 'CT'] - row[1]))

	# Calculate the SEM for technical replicate groups
	targets = set(data['Target Name'])
	mean_sem_result = {}
	for target in targets:
		mean_sem_result[target] = {}
		samples = set(data[data['Target Name'] == target]['Sample Name'])
		for sample in samples:
			target_sample_data = data[(data['Target Name'] == target) & (data['Sample Name'] == sample)]
			mean = target_sample_data['rq'].mean()
			sdt_dev = target_sample_data['rq'].std()
			std_err = target_sample_data['rq'].sem()
			mean_sem_result[target][sample] = (mean, sdt_dev, std_err)
	for i_row, row in data.iterrows():
		if data.at[i_row, 'Sample Name'] in samples and data.at[i_row, 'Sample Name'] in samples and data.at[i_row, 'Target Name'] in mean_sem_result and data.at[i_row, 'Sample Name'] in mean_sem_result[data.at[i_row, 'Target Name']]:
			data.at[i_row, 'rq'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][0]
			data.at[i_row, 'rqSD'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][1]
			data.at[i_row, 'rqSEM'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][2]

	data_output_summary = data[(data['Control'].eq(False))].groupby(['Target Name','Sample Name']).agg({'rq': [np.size, 'mean'], 'rqSD': 'mean', 'rqSEM': 'mean'})

	return data_output_summary

def process_stability(data, csample):

	# Calculate CT Mean (Endogenous Control Mean) and SSD for all Controls

	control_filter = (data['Control'].eq(True))
	data_controls_grouped = data[control_filter].groupby(['Target Name','Sample Name']).agg({'CT': [np.size, 'mean', 'std']})

	print("Endogenous Control CT Means and SSD")
	print(data_controls_grouped)

	data_controls_ct = data[control_filter].groupby(['Sample Name']).agg({'CT':'mean'})

	print("Combined Endogenous Control CT Means and SSD")
	print(data_controls_ct)

	#Create deltaCT column        
	for i, row in enumerate(data_controls_ct.itertuples(name = None), 1):
		name_filter = (data['Control'].eq(False)) & (data['Sample Name'] == row[0])
		for j in data[name_filter].index:
			data.loc[j, 'deltaCT'] = data.loc[j, 'CT'] - row[1]

	#Mark the Control Samples
	data['ControlSample'] = data['Sample Name'].apply(lambda x: True if x in csample else False)
	filter_sample = (data['ControlSample'].eq(True)) & (data['Control'].eq(False))
	data_sampled = data[filter_sample].groupby(['Target Name']).agg({('deltaCT'):'mean'})

	print("Mean Control Sample Delta CT")
	print(data_sampled)

	for i, row in enumerate(data_sampled.itertuples(name = None), 1):
		target_filter = (data['Target Name'] == row[0]) & (data['Control'].eq(False))
		for j in data[target_filter].index:
			data.loc[j, 'rq'] = np.power(2, -(data.loc[j, 'deltaCT'] - row[1]))

	# Calculate the SEM for technical replicate groups
	targets = set(data['Target Name'])
	mean_sem_result = {}
	for target in targets:
		mean_sem_result[target] = {}
		samples = set(data[data['Target Name'] == target]['Sample Name'])
		for sample in samples:
			target_sample_data = data[(data['Target Name'] == target) & (data['Sample Name'] == sample)]
			mean = target_sample_data['rq'].mean()
			sdt_dev = target_sample_data['rq'].std()
			std_err = target_sample_data['rq'].sem()
			mean_sem_result[target][sample] = (mean, sdt_dev, std_err)
	for i_row, row in data.iterrows():
		if data.at[i_row, 'Sample Name'] in samples and data.at[i_row, 'Sample Name'] in samples and data.at[i_row, 'Target Name'] in mean_sem_result and data.at[i_row, 'Sample Name'] in mean_sem_result[data.at[i_row, 'Target Name']]:
			data.at[i_row, 'rq'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][0]
			data.at[i_row, 'rqSD'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][1]
			data.at[i_row, 'rqSEM'] = mean_sem_result[data.at[i_row, 'Target Name']][data.at[i_row, 'Sample Name']][2]

	data_output_summary = data[(data['Control'].eq(False))].groupby(['Target Name','Sample Name']).agg({'rq': [np.size, 'mean'], 'rqSD': 'mean', 'rqSEM': 'mean'})

	return data_output_summary


def cleanup_outliers(d, feature, cutoff, max_outliers):
	"""Function to remove outliers based on cutoff and maximum number of outliers, 
	by removing the furthest data point in each group when the standard deviation 
	is higher than the cutoff"""

	# Calculate SSD for all sample groups
	f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN')
	d1 = d[f].groupby(['Sample Name', 'Target Name']).agg({'CT': ['std']})
	f = (d1['CT']['std'] > cutoff)
	d2 = d1[f]
	if not d2.empty:
		# Mark all outliers
		for i, row in enumerate(d2.itertuples(name = None), 1):
			f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
				& (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
			dx_idx = d[f].index
			group_size = len(dx_idx)
			min_size = round(group_size * (1 - max_outliers))
			size = group_size
			if min_size < 2:
				min_size = 2
				h.log(5, 'Minimum group size must be equal or greater than 2')
			while True:
				f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
					& (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
				dx = d[f].copy()
				dxg = d[f].groupby(['Sample Name', 'Target Name']).agg({feature: [np.size, 'std', 'mean']})
				if dxg[feature]['std'].iloc[0] <= cutoff:
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
				dx['Distance'] = (dx[feature] - dxg[feature]['mean'].iloc[0])**2
				j = dx.sort_values(by = 'Distance', ascending = False).index[0]
				d['Ignore'].loc[j] = True

	return(d[(d['Ignore'].eq(False))])

def get_cmd_args():
	"""This sets up the ArgumentParser and returns the list of arguments"""


	#Creates the Argument Parser
	parser = ArgumentParser(description = "ID Lab qPCR Analysis v" + VERSION + " " + QUALITY)

	#Adds the input file argument
	parser.add_argument('-f', '--file',
				nargs = '+',
				type = FileType('r'),
				required = True)

	#Adds the output directory
	parser.add_argument('-o', '--output',
				required = True)

	#Adds the model argument, to select between the three models
	parser.add_argument('-m', '--mod', '--model',
				nargs = '?',
				choices = ['relative', 'absolute', 'stability'],
				required = True)

	#Adds the control genes argument, taking a list of gene names
	parser.add_argument('-cg', '--cgenes', '--controlgenes',
				nargs = '+',
				required = True)

	#Adds the optional control sample argument for the stability model, taking a list of sample names
	parser.add_argument('-cs', '--csample', '--controlsamples',
				nargs = '*')

	#Adds optional outlier cutoff
	parser.add_argument('-oc', '--ocutoff',
				type = float,
				default = 0.3)

	#Adds optional max outliers
	parser.add_argument('-om', '--omax',
				type = float,
				default = 0.5)

	#Adds optional encoding 
	parser.add_argument('-e', '--encoding',
				default = 'ISO-8859-1')

	#Adds optional header size
	parser.add_argument('-hd', '--header',
				default = 47)

	return vars(parser.parse_args())

#-f C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_stability/2019-05-23_133411-controlandko-line.csv -o C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_relative/ -m stability -cg GAPDH ACTB -cs H9
#	
#-f C:/Users/eddie/Documents/idlab/test/Iva/2018-09-14GAPDHP2X3.csv C:/Users/eddie/Documents/idlab/test/Iva/2018-09-18OPRD1bACTIN.csv C:/Users/eddie/Documents/idlab/test/Iva/2018-09-20PIEZO2.csv C:/Users/eddie/Documents/idlab/test/Iva/2018-10-02RET.csv -o C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_relative/ -m absolute -cg GAPDH ACTB BACT 


if __name__ == "__main__":
	main()
