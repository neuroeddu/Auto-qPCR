import pandas as pd
import AUTOqPCR
import matplotlib.pyplot as plt
import numpy as np

# figure format: font size, bar width, error bar, facecolor
fs = 60
barwidth = 0.75
error_kw = dict(lw=4, capsize=7, capthick=4)
fc = '#F0F0F0'


def plots(dataframe, model, targets, samples):
	plots = []
	# Absolute: bar plot for each gene, a grouped bar plot by samples and a grouped bar plot by genes
	if model == 'absolute':
		plot_by_samples = plt.figure(figsize=(len(samples)*len(targets)*1.2, 30))
		counter = 0
		for item in targets:
			sample = list(dataframe.loc[item, 'NormQuant']['mean'])
			x = np.arange(len(sample))
			x2 = x * len(targets) + barwidth * counter
			plot = plt.figure(figsize=(len(sample)*1.5, 25))
			# set color, width, edgecolor, etc.
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'NormSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, label=item)
			plt.xlabel(item, fontweight='bold', fontsize=fs+10, labelpad=20)
			plt.xticks([i for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
			plt.yticks(fontsize=fs)
			plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
			plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
			# set axes width
			plt.gca().spines['bottom'].set_linewidth(5)
			plt.gca().spines['left'].set_linewidth(5)
			plt.gca().spines['top'].set_visible(False)
			plt.gca().spines['right'].set_visible(False)
			plt.gca().set_facecolor(fc)
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
		plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close(plot_by_samples)
		plots.append(plot_by_samples)

		plot_by_genes = plt.figure(figsize=(len(targets) * len(samples)*1.2, 30))
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'NormQuant']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'NormSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	# Genomic stability: grouped by chromosomes and by cell lines
	elif model == 'stability':
		# plot grouped by chromosomes
		plot_by_samples = plt.figure(figsize=(len(samples)*len(targets)*1.2, 30))
		counter = 0
		for item in targets:
			sample = list(dataframe.loc[item, 'rq']['mean'])
			x = np.arange(len(sample)) * len(targets) + barwidth * counter
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(targets) + barwidth * counter / 2 for i in range(len(samples))], samples, rotation='vertical', fontsize=fs)
		plt.yticks(fontsize=fs)
		plt.xlabel('DNA Regions', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)
		# plot grouped by chromosomes
		plot_by_chrs = plt.figure()
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'rq']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_chrs)

	# relative deltaCT and delta delta CT: plots for each gene and grouped plots by samples and genes
	else:
		plot_by_samples = plt.figure(figsize=(len(targets) * len(samples)*1.2, 30))
		counter = 0
		for item in targets:
			plot = plt.figure(figsize=(len(sample)*1.5, 25))
			sample = list(dataframe.loc[item, 'rq']['mean'])
			x = np.arange(len(sample))
			x2 = x * len(targets) + barwidth * counter
			plt.bar(x, sample, yerr=list(dataframe.loc[item, 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			plt.xlabel(item, fontsize=fs+10, fontweight='bold', labelpad=20)
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
			plt.gca().set_facecolor(fc)
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
		plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
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
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)
		# grouped by genes
		plot_by_genes = plt.figure(figsize=(len(targets) * len(samples)*1.2, 30))
		counter = 0
		for item in samples:
			target = list(dataframe.loc[(slice(None), item), 'rq']['mean'])
			x = np.arange(len(target)) * len(samples) + barwidth * counter
			plt.bar(x, target, yerr=list(dataframe.loc[(slice(None), item), 'rqSEM']['mean']), align='center',
					error_kw=error_kw, width=barwidth, edgecolor='white', label=item)
			counter += 1

		plt.xticks([i * len(samples) + barwidth * counter / 2 for i in range(len(targets))], targets, rotation='horizontal', fontsize=fs)
		plt.yticks(fontsize=fs)
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
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
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	return plots


# grouped plots with endogeneous control removed by samples and genes for absolute and relative models
def plots_wo_controls(dataframe, model, targets, samples, cgenes):
	targets = [t for t in targets if t.lower() not in cgenes.lower().split(',')]
	dataframe = dataframe.loc[targets, slice(None), :]

	plots = []

	# Absolute: a grouped bar plot by genes and a grouped bar plot by cell lines
	if model == 'absolute':
		plot_by_samples = plt.figure(figsize=(len(targets) * len(samples)*1.2, 30))
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
		plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)

		plot_by_genes = plt.figure(figsize=(len(targets) * len(samples)*1.2, 30))
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
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	elif model != 'stability':
		plot_by_samples = plt.figure(figsize=(len(targets) * len(samples)*1.2, 30))
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
		plt.xlabel('Samples', fontsize=fs+10, fontweight='bold', labelpad=20)
		if model == 'relative':
			plt.ylabel(r'Relative Quantification (RQ$_{ΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		else:
			plt.ylabel(r'Relative Quantification (RQ$_{ΔΔCT}$)', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_samples)

		# grouped by genes
		plot_by_genes = plt.figure(figsize=(len(targets) * len(samples)*1.2, 30))
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
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
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
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_genes)

	return plots


# plot by user defined groups in stats
def plot_by_groups(df, model, targets, cgenes):

	# list of groups
	groups = df['Group'].drop_duplicates(keep='first').values.tolist()

	plots = []
	if model == 'absolute':
		# remove endogeneous control genes
		targets = [t for t in targets if t.lower() not in cgenes.lower().split(',')]
		# grouped by groups on the x-axis
		counter = 0
		plot_by_group = plt.figure(figsize=(len(groups)*len(targets)*2, 20))
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
		plt.xlabel('Groups', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_group)
		# grouped by targets on the x-axis
		counter = 0
		plot_by_target = plt.figure(figsize=(len(groups)*len(targets)*2, 20))
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
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Normalized Expression', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_target)
	elif model == 'stability':
		# grouped by groups on the x-axis
		counter = 0
		plot_by_group = plt.figure(figsize=(len(groups)*len(targets)*2, 2))
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
		plt.xlabel('Groups', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_group)
		# grouped by groups on the x-axis
		counter = 0
		plot_by_target = plt.figure(figsize=(len(groups)*len(targets)*2, 20))
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
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.ylabel('Copy Number per Chromosome', fontsize=fs+10, fontweight='bold', labelpad=20)
		plt.legend(fontsize=fs, loc='upper left', bbox_to_anchor=(1, 1))
		# set axes width
		plt.gca().spines['bottom'].set_linewidth(5)
		plt.gca().spines['left'].set_linewidth(5)
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_target)
	else:
		# remove endogeneous control genes
		targets = [t for t in targets if t.lower() not in cgenes.lower().split(',')]
		# grouped by groups on the x-axis
		counter = 0
		plot_by_group = plt.figure(figsize=(len(groups)*len(targets)*2, 20))
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
		plt.xlabel('Groups', fontsize=fs+10, fontweight='bold', labelpad=20)
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
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_group)
		# grouped by groups on the x-axis
		counter = 0
		plot_by_target = plt.figure(figsize=(len(groups)*len(targets)*2, 20))
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
		plt.xlabel('Targets', fontsize=fs+10, fontweight='bold', labelpad=20)
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
		plt.gca().set_facecolor(fc)
		plt.gca().tick_params(width=5)
		plt.tight_layout()
		plt.close()
		plots.append(plot_by_target)

	return plots
