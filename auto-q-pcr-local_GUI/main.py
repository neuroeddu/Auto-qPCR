from flask import Flask, request, make_response, render_template, flash, redirect, url_for
import io
import pandas as pd
import AUTOqPCR
import plot
import statistics
import re
import datetime
from zipfile import ZipFile
import psutil
import os
import logging


app: Flask = Flask(__name__)
app.secret_key = b'_5#y2L"F4Q8z\n\xec]/'


@app.route('/')
def form():
	return render_template('layout.html')


@app.route('/hello')
def hello():
	return render_template('hello.html')


@app.route('/download', methods=["POST"])
def transform_view():
	# make log file
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)
	log_stream = io.StringIO()
	stream_handler = logging.StreamHandler(log_stream)
	logger.addHandler(stream_handler)

	# get current machine time
	now = datetime.datetime.now()
	date_string = now.strftime("%m-%d-%Y")

	try:
		logger.info('Started')
		files = request.files.getlist('file[]')
		if len(files) == 0:
			return "No file"
		# Creates empty data frame
		data = pd.DataFrame()

		for item in files:
			stream = io.StringIO(item.stream.read().decode("utf-8"), newline="")
			i = 0
			while True:
				h = stream.readline()
				if not h:
					i = -1
					break
				# if not h.replace(",",""):
				# i += 1
				if h.upper().startswith("WELL"):
					break
				i += 1
			if i == -1:
				print("No header found in '{}'".format(item))
				continue
			else:
				print("Header found at {} in '{}'".format(i, item))
			# print(i)

			stream.seek(0)
			filedata = pd.read_csv(stream,
								   skip_blank_lines=True,
								   skipinitialspace=True,
								   engine='python',
								   encoding="utf-8",
								   header=i)

			# print(filedata)
			filedata['filename'] = item.filename
			filedata.rename(columns=rx_rename, inplace=True)
			data = data.append(filedata, ignore_index=True, sort=True)

		# stream.seek(0)
		logger.info('Files upload complete.')

		genes = request.form['genes']
		logger.info('Gene names if they are included in file names: ' + genes)

		if genes != '':
			genes = [genes.strip() for genes in genes.split(',')]
			data['Target Name'] = data['filename'].str.extract(re.compile('(' + '|'.join(genes) + ')', re.IGNORECASE),
															   expand=False).fillna('')
		model = request.form['option']
		quencher = request.form['quencher']
		task = request.form['task']
		cgenes = request.form['cgenes']
		cutoff = request.form.get('cutoff', type=float)
		max_outliers = request.form.get('max_outliers', type=float)
		target_sorter = request.form['target_sorter']
		sample_sorter = request.form['sample_sorter']
		csample = request.form['csample']
		colnames = request.form['colnames']
		qty = request.form.get('quantity', type=int)
		gcol = request.form['gcol']
		glist = request.form['glist']
		rm = request.form['option2']
		nd = request.form['option4']

		logger.info('Model: ' + model + '\nEndogenous control genes: ' + cgenes + '\nCut-off: ' + str(cutoff) + \
					'\nMaximum Outliers: ' + str(max_outliers) + '\nTarget Order: ' + target_sorter + '\nSample Order: ' + \
					sample_sorter + '\nControl Sample: ' + csample + '\nAdditional column names: ' + colnames + \
					'\nNumber of groups: ' + str(qty) + '\nGroup column name: ' + gcol + '\nGroup name: ' + glist + \
					'\nRepeated measures: ' + rm + '\n' + 'Normal distribution: ' + nd)

		clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(data, model, quencher, task, cgenes, cutoff,
																		   max_outliers,
																		   target_sorter, sample_sorter, csample, colnames)
		# making summary data csv
		output = summary_data.to_csv()
		output_w_group = summary_data_w_group.to_csv()
		clean_output = clean_data.to_csv()

		logger.info('Clean data and summary data are created')

		plots = plot.plots(summary_data, model, targets, samples)
		plots2 = plot.plots_wo_controls(summary_data, model, targets, samples, cgenes)

		logger.info('Plots of the summary data are created.')

		# making stats csv
		if qty is not None:
			if request.form['option3'] != 'False':
				if gcol.lower() in clean_data.columns.str.lower():
					col = gcol
				clean_data['Group'] = clean_data[col]
			else:
				groups = glist.split(',')
				clean_data = statistics.add_groups(clean_data, groups)

			stats_dfs, posthoc_dfs = statistics.stats(model, qty, clean_data, targets, rm, nd)
			stats_output = stats_dfs.to_csv(index=False)
			posthoc_output = posthoc_dfs.to_csv(index=False)

			logger.info('Statistics output data are created.')

			group_plot = plot.plot_by_groups(clean_data, model, targets, cgenes)

			logger.info('Plots of statistics output are created.')

		outfile = io.BytesIO()
		with ZipFile(outfile, 'w') as myzip:
			if qty is not None:
				if qty == 2:
					if nd == 'True':
						myzip.writestr('ttest_result.csv', stats_output)
					else:
						if rm == 'True':
							myzip.writestr('MannWhitneyUTest_result.csv', stats_output)
						else:
							myzip.writestr('WilcoxonTest_result.csv', stats_output)
				else:
					if nd == 'True':
						myzip.writestr('ANOVA_result.csv', stats_output)
						myzip.writestr('Posthoc_result.csv', posthoc_output)
					else:
						if rm == 'True':
							myzip.writestr('Friedman_result.csv', stats_output)
						else:
							myzip.writestr('KruskalWallisTest_result.csv', stats_output)
						myzip.writestr('Posthoc_result.csv', posthoc_output)

				buf = io.BytesIO()
				group_plot[0].savefig(buf)
				image_name = 'Plot_by_groups.png'
				myzip.writestr(image_name, buf.getvalue())
				buf.close()
				buf = io.BytesIO()
				group_plot[1].savefig(buf)
				image_name2 = 'Plot_by_targets.png'
				myzip.writestr(image_name2, buf.getvalue())
				buf.close()
			myzip.writestr('clean_data.csv', clean_output)
			myzip.writestr('summary_data.csv', output)
			myzip.writestr('summary_data_w_groups.csv', output_w_group)
			myzip.writestr('log.txt', log_stream.getvalue())
			log_stream.flush()
			# individual plots
			for i in range(len(plots) - 2):
				buf = io.BytesIO()
				plots[i].savefig(buf)
				if model != 'instability':
					image_name = targets[i] + '.png'
					myzip.writestr(image_name, buf.getvalue())
			# grouped plots by sample and by genes
			buf = io.BytesIO()
			plots[len(plots) - 2].savefig(buf)
			image_name = 'Sample_Groups.png'
			myzip.writestr(image_name, buf.getvalue())
			buf.close()
			buf = io.BytesIO()
			plots[len(plots) - 1].savefig(buf)
			image_name = 'All_Targets.png'
			myzip.writestr(image_name, buf.getvalue())
			buf.close()
			if model != 'instability':
				# plots with endogeneous controls removed
				buf = io.BytesIO()
				plots2[0].savefig(buf)
				image_name = 'Sample_Groups (without endogenous controls).png'
				myzip.writestr(image_name, buf.getvalue())
				buf.close()
				buf = io.BytesIO()
				plots2[1].savefig(buf)
				image_name = 'All_Targets (without endogenous controls).png'
				myzip.writestr(image_name, buf.getvalue())
				buf.close()
				myzip.close()

		response = make_response(outfile.getvalue())
		response.headers['Content-Type'] = 'application/zip'
		response.headers['Content-Disposition'] = 'attachment; filename=outputs_' + model + '_' + date_string + '.zip'
		outfile.close()
		# alert
		flash('Your data has been processed successfully!', 'success')
	except Exception as e:
		logger.error('Error occurred: ' + str(e))
		response = make_response(log_stream.getvalue())
		response.headers['Content-Type'] = 'text/plain'
		response.headers['Content-Disposition'] = 'attachment; filename=log_'+date_string+'.txt'
		log_stream.flush()
		# alert
		flash('Sorry, something went wrong. Please check log.txt file.', 'danger')

	# get CPU and memory for the process
	# myProcess = psutil.Process(os.getpid())
	# print('CPU percent: ' + str(myProcess.cpu_percent()))
	# print('Memory info: ' + str(myProcess.memory_info()[0]/2.**30))

	return response


def rx_rename(col_name):
	# compiles regular expressions
	rct = re.compile('((?<![\w _])[(]*c[()ycle ]*t[)hreshold]*(?![\w\W]))', re.IGNORECASE)
	# rquant = re.compile('((?<![\w _])[(]*quant[ity)]*\Z(?! sd)(?! mean))|(?<![\w _])(ng)', re.IGNORECASE)
	rquant = re.compile('((?<![\w _])[(]*quant[ity) ]*\Z(?! sd)(?! mean))', re.IGNORECASE)
	rsamp = re.compile('(samp)+|(chrom)+|(desc)+', re.IGNORECASE)
	rtarg = re.compile('(targ)+|(gene)+', re.IGNORECASE)
	rdye = re.compile('(repor)+|(dye)+|(fluor)+', re.IGNORECASE)
	rtask = re.compile('(task)+|(role)+|(content)+', re.IGNORECASE)

	if re.match(rct, col_name):
		return str('CT')
	if re.match(rquant, col_name):
		return str('Quantity')
	if re.match(rsamp, col_name):
		return str('Sample Name')
	if re.match(rtarg, col_name):
		return str('Target Name')
	if re.match(rdye, col_name):
		return str('Reporter')
	if re.match(rtask, col_name):
		return str('Task')
	else:
		return col_name


if __name__ == '__main__':
	app.debug = True
	app.run(host='0.0.0.0', port=5000)
