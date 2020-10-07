import pandas as pd
from application import AUTOqPCR
import matplotlib.pyplot as plt
import numpy as np

# figure format: font size, bar width, error bar, facecolor
fs = 80
barwidth = 0.75
error_kw = dict(lw=4, capsize=7, capthick=4)
# fc = '#F0F0F0'


def plots(dataframe, model, targets, samples):

	if len(samples)*len(targets) < 5:
		figsize = (20, 25)
	elif len(samples)*len(targets) < 10:
		figsize = (len(samples)*len(targets)*5, 25)
	elif len(samples)*len(targets) < 25:
		figsize = (len(samples)*len(targets)*3.2, 25)
	else:
		figsize = (len(samples)*len(targets)*1.6, 25)

	if len(samples) < 20:
		sfigsize = (len(samples) * 5, 27)
	else:
		sfigsize = (len(samples) * 2, 27)

	# number of columns in the legend
	ncol_s = len(samples) // 10 + 1
	ncol_t = len(targets) // 10 + 1

	plots = []
	# Absolute: bar plot for each gene, a grouped bar plot by samples and a grouped bar plot by genes
	if model == 'absolute':
		plot_by_samples = plt.figure(figsize=figsize)
		counter = 0
		for item in targets:
			sample = list(dataframe.loc[item, 'NormQuant']['mean'])
			x = np.arange(len(samples))
			x2 = x * len(targets) + barwidth * counter
			plot = plt.figure(figsize=sfigsize)
			# set color, width, edgecolor, etc.
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'NormSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, label=item)
			# plt.xlabel(item, fontweight='bold', fontsize=fs+10, labelpad=20)
			plt.xticks([i for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
			plt.yticks(fontsize=fs)
			plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close(plot)
			plots.append(plot)
			# grouped plot
			plt.bar(x2, sample, yerr=list(dataframe.loc[item, 'NormSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1
		plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)

		plot_by_genes = plt.figure(figsize=figsize)
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'NormQuant']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'NormSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_s)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	# Genomic stability: grouped by chromosomes and by cell lines
	elif model == 'stability':
		# plot grouped by DNA regions
		plot_by_samples = plt.figure(figsize=figsize)
		counter = 0
		for item in targets:
			sample = list(dataframe.loc[item, 'rq']['mean'])
			x = np.arange(len(sample)) * len(targets) + barwidth * counter
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('DNA Regions', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)
		# plot grouped by chromosomes
		plot_by_chrs = plt.figure(figsize=figsize)
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'rq']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_s)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_chrs)

	# relative deltaCT and delta delta CT: plots for each gene and grouped plots by samples and genes
	else:
		plot_by_samples = plt.figure(figsize=figsize)
		counter = 0
		for item in targets:
			plot = plt.figure(figsize=sfigsize)
			sample = list(dataframe.loc[item, 'rq']['mean'])
			x = np.arange(len(samples))
			x2 = x * len(targets) + barwidth * counter
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			# plt.xlabel(item, fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.xticks([i for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
			plt.yticks(fontsize=fs)

			if model == 'relative_dCT':
				plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
			else:
				plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close(plot)
			plots.append(plot)
			# grouped plot by samples
			plt.bar(x2, sample, yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1
		plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
		if model == 'relative_dCT':
			plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		else:
			plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		# plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)
		# grouped by genes
		plot_by_genes = plt.figure(figsize=figsize)
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'rq']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		if model == 'relative_dCT':
			plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		else:
			plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_s)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	return plots


# grouped plots with endogeneous control removed by samples and genes for absolute and relative models
def plots_wo_controls(dataframe, model, targets, samples, cgenes):
	targets = [t for t in targets if t.lower() not in cgenes.lower().split(',')]
	dataframe = dataframe.loc[targets, slice(None), :]

	# number of columns in the legend
	ncol_s = len(samples) // 10 + 1
	ncol_t = len(targets) // 10 + 1

	#set figure size
	if len(samples) * len(targets) < 5:
		figsize = (20, 25)
	elif len(samples) * len(targets) < 10:
		figsize = (len(samples) * len(targets) * 5, 25)
	elif len(targets) * len(samples)*1.6 < 40:
		figsize = (len(targets) * len(samples)*3.2, 25)
	else:
		figsize = (len(targets) * len(samples)*1.6, 25)

	plots = []

	# Absolute: a grouped bar plot by genes and a grouped bar plot by cell lines
	if model == 'absolute':
		plot_by_samples = plt.figure(figsize=figsize)
		counter = 0
		for item in targets:
			sample = list(dataframe.loc[item, 'NormQuant']['mean'])
			x = np.arange(len(sample)) * len(targets) + barwidth * counter
			# grouped plot
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'NormSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white',
					label=item)
			counter += 1

		plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)

		plot_by_genes = plt.figure(figsize=figsize)
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'NormQuant']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'NormSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white',
					label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_s)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	elif model != 'stability':
		plot_by_samples = plt.figure(figsize=figsize)
		counter = 0
		for item in targets:
			sample = list(dataframe.loc[item, 'rq']['mean'])
			x = np.arange(len(sample)) * len(targets) + barwidth * counter
			# grouped plot by samples
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white',
					label=item)
			counter += 1
		plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
		if model == 'relative':
			plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		else:
			plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)

		# grouped by genes
		plot_by_genes = plt.figure(figsize=figsize)
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'rq']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white',
					label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		# plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		if model == 'relative_dCT':
			plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		else:
			plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_s)
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		# plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	return plots


# plot by user defined groups in stats
def plot_by_groups(df, model, targets, cgenes, tw):
	if tw == 'False':
		# list of groups
		groups = df['Group'].drop_duplicates(keep='first').values.tolist()

		# set figure size:
		if len(groups) * len(targets) * 2.4 < 10:
			figsize = (len(groups) * len(targets) * 4.5, 20)
		elif len(groups) * len(targets) * 2.4 < 20:
			figsize = (len(groups) * len(targets) * 3.4, 20)
		else:
			figsize = (len(groups) * len(targets) * 2.4, 20)

		# number of columns in the legend
		ncol_t = len(targets) // 10 + 1
		ncol_g = len(groups) // 10 + 1

		plots = []
		if model == 'absolute':
			# remove endogeneous control genes
			targets = [t for t in targets if t.lower() not in cgenes.lower().split(',')]
			# grouped by groups on the x-axis
			counter = 0
			plot_by_group = plt.figure(figsize=figsize)
			for t in targets:
				y = []
				st_err = []
				x = np.arange(len(groups)) * len(targets) + barwidth * counter
				for g in groups:
					sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
					y.append(sample['NormMean'].mean())
					st_err.append(sample['NormMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=t)

				counter += 1
			plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(groups))], groups, rotation='vertical', fontsize=fs)
			plt.yticks(fontsize=fs)
			# plt.xlabel('Groups', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
			plots.append(plot_by_group)
			# grouped by targets on the x-axis
			counter = 0
			plot_by_target = plt.figure(figsize=figsize)
			for g in groups:
				y = []
				st_err = []
				x = np.arange(len(targets)) * len(groups) + barwidth * counter
				for t in targets:
					sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
					y.append(sample['NormMean'].mean())
					st_err.append(sample['NormMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=g)
				counter += 1
			plt.xticks([i * len(groups) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
			plt.yticks(fontsize=fs)
			# plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_g)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
			plots.append(plot_by_target)
		elif model == 'stability':
			# grouped by groups on the x-axis
			counter = 0
			plot_by_group = plt.figure(figsize=figsize)
			for t in targets:
				y = []
				st_err = []
				x = np.arange(len(groups)) * len(targets) + barwidth * counter
				for g in groups:
					sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
					y.append(sample['rqMean'].mean())
					st_err.append(sample['rqMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=t)

				counter += 1
			plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(groups))], groups, rotation='horizontal', fontsize=fs)
			plt.yticks(fontsize=fs)
			# plt.xlabel('Groups', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
			plots.append(plot_by_group)
			# grouped by groups on the x-axis
			counter = 0
			plot_by_target = plt.figure(figsize=figsize)
			for g in groups:
				y = []
				st_err = []
				x = np.arange(len(targets)) * len(groups) + barwidth * counter
				for t in targets:
					sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
					y.append(sample['rqMean'].mean())
					st_err.append(sample['rqMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=g)

				counter += 1
			plt.xticks([i * len(groups) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
			plt.yticks(fontsize=fs)
			# plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_g)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
			plots.append(plot_by_target)
		else:
			# remove endogeneous control genes
			targets = [t for t in targets if t.lower() not in cgenes.lower().split(',')]
			# grouped by groups on the x-axis
			counter = 0
			plot_by_group = plt.figure(figsize=figsize)
			for t in targets:
				y = []
				st_err = []
				x = np.arange(len(groups)) * len(targets) + barwidth * counter
				for g in groups:
					sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
					y.append(sample['rqMean'].mean())
					st_err.append(sample['rqMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=t)

				counter += 1
			plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(groups))], groups, rotation='horizontal', fontsize=fs)
			plt.yticks(fontsize=fs)
			# plt.xlabel('Groups', fontsize=fs+10, fontweight='bold', labelpad=20)
			if model == 'relative_dCT':
				plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
			else:
				plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_t)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
			plots.append(plot_by_group)
			# grouped by groups on the x-axis
			counter = 0
			plot_by_target = plt.figure(figsize=figsize)
			for g in groups:
				y = []
				st_err = []
				x = np.arange(len(targets)) * len(groups) + barwidth * counter
				for t in targets:
					sample = df.loc[(df['Target Name'] == t) & (df['Group'] == g)]
					y.append(sample['rqMean'].mean())
					st_err.append(sample['rqMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=g)

				counter += 1
				plt.xticks([i * len(groups) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
				plt.yticks(fontsize=fs)
			if model == 'relative_dCT':
				plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
			else:
				plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_g)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
			plots.append(plot_by_target)
	else:
		group1 = df['Group1'].drop_duplicates(keep='first').values.tolist()
		group2 = df['Group2'].drop_duplicates(keep='first').values.tolist()

		# number of columns in the legend
		ncol_g1 = len(group1) // 10 + 1

		# set figure size:
		if len(group1) * len(group2) * 2.4 < 10:
			figsize = (len(group1) * len(group2) * 4.5, 20)
		elif len(group1) * len(group2) * 2.4 < 20:
			figsize = (len(group1) * len(group2) * 3.4, 20)
		else:
			figsize = (len(group1) * len(group2) * 2.4, 20)

		plots = []

		if model == 'absolute':
			# grouped by groups on the x-axis
			counter = 0
			plot = plt.figure(figsize=figsize)
			for g1 in group1:
				y = []
				st_err = []
				x = np.arange(len(group2)) * len(group1) + barwidth * counter
				for g2 in group2:
					sample = df.loc[(df['Group1'] == g1) & (df['Group2'] == g2)]
					y.append(sample['NormMean'].mean())
					st_err.append(sample['NormMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=g1)

				counter += 1
			plt.xticks([i * len(group1) + barwidth * counter / 2 for i in range(len(group2))], group2, rotation='vertical', fontsize=fs)
			plt.yticks(fontsize=fs)
			plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_g1)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
		elif model == 'stability':
			# grouped by groups on the x-axis
			counter = 0
			plot = plt.figure(figsize=figsize)
			for g1 in group1:
				y = []
				st_err = []
				x = np.arange(len(group2)) * len(group1) + barwidth * counter
				for g2 in group2:
					sample = df.loc[(df['Group1'] == g1) & (df['Group2'] == g2)]
					y.append(sample['rqMean'].mean())
					st_err.append(sample['rqMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=g1)

				counter += 1
			plt.xticks([i * len(group1) + barwidth * counter / 2 for i in range(len(group2))], group2, rotation='vertical', fontsize=fs)
			plt.yticks(fontsize=fs)
			plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_g1)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
		else:
			# grouped by groups on the x-axis
			counter = 0
			plot = plt.figure(figsize=figsize)
			for g1 in group1:
				y = []
				st_err = []
				x = np.arange(len(group2)) * len(group1) + barwidth * counter
				for g2 in group2:
					sample = df.loc[(df['Group1'] == g1) & (df['Group2'] == g2)]
					y.append(sample['rqMean'].mean())
					st_err.append(sample['rqMean'].sem())
				plt.bar(x, y, yerr=st_err, error_kw=error_kw, align='center',
						width=barwidth, edgecolor='white', label=g1)

				counter += 1
			plt.xticks([i * len(group1) + barwidth * counter / 2 for i in range(len(group2))], group2,
					   rotation='vertical', fontsize=fs)
			plt.yticks(fontsize=fs)
			if model == 'relative_dCT':
				plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs + 10, fontweight='bold', labelpad=20)
			else:
				plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs + 10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1), ncol=ncol_g1)
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			# plt.gca().set_facecolor(fc)
			plt.gca().tick_params(width=5)
			plt.tight_layout()
			plt.close()
		plots.append(plot)

	return plots
