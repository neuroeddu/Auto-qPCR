from flask import Flask, request, make_response, render_template
import io
import pandas as pd
import AUTOqPCR
import plot
import statistics
import re
from zipfile import ZipFile
import datetime

app: Flask = Flask(__name__)


@app.route('/')
def form():
	return render_template('layout.html')


@app.route('/download', methods=["POST"])
def transform_view():

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
			#if not h.replace(",",""):
			#i += 1
			if h.upper().startswith("WELL"):
				break
			i += 1
		if i == -1:
			print("No header found in '{}'".format(item))
			continue
		else:
			print("Header found at {} in '{}'".format(i, item))
		#print(i)

		stream.seek(0)
		filedata = pd.read_csv(stream ,
							   skip_blank_lines=True ,
							   skipinitialspace=True ,
							   engine='python' ,
							   encoding="utf-8" ,
							   header= i)

		# print(filedata)	
		data = data.append(filedata , ignore_index=True , sort=True)
		data['filename'] = item.filename
		#stream.seek(0)

	model = request.form['option']
	cgenes = request.form['cgenes']
	cutoff = request.form.get('cutoff' , type=float)
	max_outliers = request.form.get('max_outliers', type=float)
	target_sorter = request.form['target_sorter']
	sample_sorter = request.form['sample_sorter']
	csample = request.form['csample']
	qty = request.form.get('quantity', type=int)
	rm = request.form['option2']
	nd = request.form['option5']
	posthoc = request.form['option3']

	data1, summary_data, targets, samples = AUTOqPCR.process_data(data, model, cgenes, cutoff, max_outliers, target_sorter, sample_sorter, csample)

	plots = plot.plots(summary_data , model , targets , samples, cgenes)
	
	# making stats csv
	if qty is not None:
		if request.form['option4'] != 'False':
			gcol = data[request.form['gcol']]
			data1['Group'] = gcol
		else:
			groups = request.form['glist'].split(',')
			data1 = statistics.add_groups(data1, groups)
		# print(data1)
		group_plot = plot.plot_by_groups(data1, model, targets, groups)

		stats_dfs, posthoc_dfs = statistics.stats(model, qty, data1, targets, rm, nd, posthoc)
		stats_output = stats_dfs.to_csv(index=False)
		posthoc_output = posthoc_dfs.to_csv(index=False)

	# making summary data csv
	output = summary_data.to_csv()
	clean_output = data1.to_csv()

	# remove cgenes from targets
	targets2 = [g for g in targets if g.lower() not in cgenes.lower().split(',')]

	# get current machine time
	now = datetime.datetime.now()
	date_string = now.strftime("%m-%d-%Y")

	#get name for output
	model_map = {'absolute':'absolute', 'relative':'relative_dCT', 'stability':'relative_ddCT', 'stability2':'genomic_stability'}
	model_name = model_map[model]

	outfile = io.BytesIO()
	with ZipFile(outfile, 'w') as myzip:
		if qty is not None:
			if qty == 2:
				if nd == 'False':
					myzip.writestr('ttest_result.csv', stats_output)
				else:
					myzip.writestr('MannWhitneyUTest_result.csv' , stats_output)
			else:
				if nd == 'False':
					myzip.writestr('anova_result.csv' , stats_output)
				else:
					myzip.writestr('KruskalWallisTest_result.csv' , stats_output)
				myzip.writestr(posthoc+'_result.csv' , posthoc_output)

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
		for i in range(len(plots)-2):
			buf = io.BytesIO()
			plots[i].savefig(buf)
			if model != 'stability2':
				image_name = targets2[i]+'.png'
			myzip.writestr(image_name, buf.getvalue())
		buf = io.BytesIO()
		plots[len(plots) - 2].savefig(buf)
		image_name = 'All_Genes.png'
		myzip.writestr(image_name , buf.getvalue())
		buf.close()
		buf = io.BytesIO()
		plots[len(plots) - 1].savefig(buf)
		image_name = 'Sample_Groups.png'
		myzip.writestr(image_name , buf.getvalue())
		buf.close()
		myzip.close()

	response = make_response(outfile.getvalue())
	response.headers['Content-Type'] = 'application/actet-stream'
	response.headers['Content-Disposition'] = 'attachment; filename=outputs_'+model_name+'_' + date_string + '.zip'
	outfile.close()

	return response


if __name__=='__main__':
	app.debug = True
	app.run(host = '0.0.0.0', port=5000)

