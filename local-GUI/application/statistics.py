import pandas
import pingouin as pg
from pingouin import pairwise_ttests, multicomp, ttest
from scipy.stats import mannwhitneyu, wilcoxon
import re


def stats(model, quantity, data, targets, tw, rm, nd):
	if model == 'absolute':
		data = data.drop(['NormQuant'], axis=1)
		data['NormMean'] = data['NormMean'].astype(float)
		mean = 'NormMean'
	else:
		data = data.drop(['rq'], axis=1)
		data['rqMean'] = data['rqMean'].astype(float)
		mean = 'rqMean'

	# prepare data from intermediate dataframe
	data = data[data['Outliers'].eq(False)]
	data = data.drop_duplicates(keep='first')

	# t-test and anova for normally distributed data
	if nd == 'True':
		if quantity == 2:
			# T-Test between 2 groups
			stats_dfs = pandas.DataFrame()
			posthoc_dfs = pandas.DataFrame()
			group = data['Group'].dropna()
			group = group.drop_duplicates(keep='first').values.tolist()
			for item in targets:
				df = data[data['Target Name'].eq(item)]
				group1 = df[df['Group'].eq(group[0])][mean]
				group2 = df[df['Group'].eq(group[1])][mean]
				t_test = ttest(group1, group2, paired=bool(rm))

				if rm == 'True':
					t_test['paired'] = 'TRUE'
				else:
					t_test['paired'] = 'FALSE'
				t_test['Target Name'] = item
				if stats_dfs is None:
					stats_dfs = t_test
				else:
					stats_dfs = stats_dfs.append(t_test, ignore_index=True)
			# reformat output table
			stats_dfs = stats_dfs.rename(columns={'cohen-d': 'effect size', 'BF10': 'Bayes factor', 'dof': 'DF'})
			cols = ['Target Name', 'DF', 'T', 'tail', 'paired', 'p-val', 'effect size', 'power', 'Bayes factor']
			stats_dfs = stats_dfs.reindex(columns=cols)
		elif quantity >= 3:
			# ANOVA test
			stats_dfs = pandas.DataFrame()
			posthoc_dfs = pandas.DataFrame()
			# tukey_dfs = pandas.DataFrame()
			pvals = []
			for item in targets:
				if rm == 'True':
					# one-way
					if tw == 'False':
						# repeated measure anova
						aov = pg.rm_anova(dv=mean, data=data[data['Target Name'].eq(item)], within='Group', subject='Sample Name', detailed=True)
						pvals.append(aov['p-unc'][0])
						aov = aov.drop([1])
						aov['measures'] = ['dependent']
						aov['Target Name'] = item
					# two-way
					else:
						aov = pg.rm_anova(dv=mean, data=data[data['Target Name'].eq(item)], within=['Group1', 'Group2'],
										  subject='Sample Name', detailed=True)
						reject_tw, pval_corr_tw = pg.multicomp(list(aov['p-unc']), alpha=0.05, method='bonf')
						aov['p-value corrected'] = pval_corr_tw
						aov['measures'] = ['dependent']*3
						aov['Target Name'] = [item]*3
					aov.drop(['eps'], axis=1)
					ph = pairwise_ttests(data=data[data['Target Name'].eq(item)], dv=mean, within='Group',
											 subject='Sample Name', padjust='fdr_bh')
					ph['Target Name'] = item
					ph['Test'] = 'T-Test'
				else:
					# one-way
					if tw == 'False':
						aov = pg.anova(dv=mean, between='Group', data=data[data['Target Name'].eq(item)], detailed=True)
						pvals.append(aov['p-unc'][0])
						aov = aov.drop([1])
						aov['measures'] = ['independent']
						aov['Target Name'] = item
						ph = pairwise_ttests(data=data[data['Target Name'].eq(item)], dv=mean, between='Group',
											 padjust='fdr_bh')
						ph['Test'] = 'T-Test'
					# two-way
					else:
						aov = pg.anova(dv=mean, between=['Group1', 'Group2'], data=data[data['Target Name'].eq(item)], detailed=False)
						aov = aov.drop([3])
						reject_tw, pval_corr_tw = pg.multicomp(list(aov['p-unc']), alpha=0.05, method='bonf')
						aov['p-value corrected'] = pval_corr_tw
						aov['measures'] = ['independent']*3
						aov['Target Name'] = [item]*3
						ph = pairwise_ttests(data=data[data['Target Name'].eq(item)], dv=mean, between=['Group1', 'Group2'], padjust='fdr_bh')
						ph['Test'] = 'T-Test'
					ph['Target Name'] = item
				if stats_dfs is None:
					stats_dfs = aov
				else:
					stats_dfs = stats_dfs.append(aov, ignore_index=True)
				if posthoc_dfs is None:
					posthoc_dfs = ph
				else:
					posthoc_dfs = posthoc_dfs.append(ph, ignore_index=True)

			reject, pvals_corr = pg.multicomp(pvals, alpha=0.05, method='bonf')

			# reformat output tables
			stats_dfs = stats_dfs.rename(columns={'p-unc': 'p-value', 'np2': 'effect size'})
			if tw == 'False':
				stats_dfs['p-value corrected'] = pvals_corr
				stats_dfs['distribution'] = ['parametric'] * len(targets)
				stats_dfs['test'] = ['ANOVA'] * len(targets)
				stats_dfs['statistic'] = ['NA'] * len(targets)
			else:
				stats_dfs['distribution'] = ['parametric'] * (len(targets)*3)
				stats_dfs['test'] = ['ANOVA'] * (len(targets)*3)
				stats_dfs['statistic'] = ['NA'] * (len(targets)*3)
			cols = ['Target Name', 'Source', 'DF', 'F', 'MS', 'SS', 'p-value', 'p-value corrected', 'measures', 'distribution', 'test',
					'statistic', 'effect size']
			stats_dfs = stats_dfs.reindex(columns=cols)

			posthoc_dfs = posthoc_dfs.drop(['Contrast', 'T'], axis=1)
			posthoc_dfs = posthoc_dfs.rename(columns={'hedges': 'effect size', 'p-corr': 'p-value corrected', 'p-unc': 'p-value',
													  'p-adjust': 'correction method', 'BF10': 'Bayes factor', 'dof': 'DF'})
			cols2 = ['Target Name', 'A', 'B', 'DF', 'p-value corrected', 'p-value', 'correction method', 'Paired',
					 'Parametric', 'Test', 'effect size', 'Bayes factor']
			posthoc_dfs = posthoc_dfs.reindex(columns=cols2)

	# nonparametric tests for not normally distributed data
	else:
		if quantity == 2:
			stats_dfs = pandas.DataFrame()
			posthoc_dfs = pandas.DataFrame()
			group = data['Group'].dropna()
			group = group.drop_duplicates(keep='first').values.tolist()
			for item in targets:
				df = data[data['Target Name'].eq(item)]
				group1 = df[df['Group'].eq(group[0])][mean]
				group2 = df[df['Group'].eq(group[1])][mean]
				if rm == 'True':
					# Mann-Whitney U test
					test = mannwhitneyu(group1, group2)
					test = pandas.DataFrame({'Target Name': item, 'pvalue': test.pvalue, 'statistic': test.statistic}, index=[0])
				else:
					# Wilcoxon
					test = wilcoxon(group1, group2)
					test = pandas.DataFrame({'Target Name': item, 'pvalue': test.pvalue, 'statistic': test.statistic}, index=[0])
				if stats_dfs is None:
					stats_dfs = test
				else:
					stats_dfs = stats_dfs.append(test, ignore_index=True)

		elif quantity >= 3:
			stats_dfs = pandas.DataFrame()
			posthoc_dfs = pandas.DataFrame()

			pvals = []
			for item in targets:
				if rm == 'True':
					# friedman test for repeated measurements
					df = pg.friedman(dv=mean, within='Group', subject='Sample Name', data=data[data['Target Name'].eq(item)])
					pvals.append(df['p-unc'][0])
					df['test'] = ['Friedman Q']
					df['measures'] = ['dependent']
					df = df.rename(columns={'Q': 'statistic'})
					df['Target Name'] = item
					df['DF'] = 'NA'
					ph = pairwise_ttests(data=data[data['Target Name'].eq(item)], dv=mean, within='Group',
										 subject='Sample Name', padjust='fdr_bh', parametric=False)
					ph['Target Name'] = item
					ph['DF'] = 'NA'
					ph['Bayes factor'] = 'NA'
					ph['Test'] = 'Wilcoxon'
				else:
					# Kruskal-Wallis H test
					df = pg.kruskal(dv=mean, between='Group', data=data[data['Target Name'].eq(item)])
					pvals.append(df['p-unc'][0])
					df['test'] = ['Kruskal-Wallis H']
					df['measures'] = ['independent']
					df = df.rename(columns={'H': 'statistic'})
					df['Target Name'] = item
					df['DF'] = 'NA'
					ph = pairwise_ttests(data=data[data['Target Name'].eq(item)], dv=mean, between='Group',
										 padjust='fdr_bh', parametric=False)
					ph['Target Name'] = item
					ph['DF'] = 'NA'
					ph['Bayes factor'] = 'NA'
					ph['Test'] = 'Mann-Whitney U'
				if stats_dfs is None:
					stats_dfs = df
				else:
					stats_dfs = stats_dfs.append(df, ignore_index=True)
				if posthoc_dfs is None:
					posthoc_dfs = ph
				else:
					posthoc_dfs = posthoc_dfs.append(ph, ignore_index=True)

			reject, pvals_corr = pg.multicomp(pvals, alpha=0.05, method='bonf')
			# reformat output tables
			stats_dfs = stats_dfs.rename(columns={'dof': 'DF', 'p-unc': 'p-value'})
			stats_dfs['p-value corrected'] = pvals_corr
			stats_dfs['distribution'] = ['non-parametric'] * len(targets)
			stats_dfs['MS'] = ['NA'] * len(targets)
			stats_dfs['SS'] = ['NA'] * len(targets)
			stats_dfs['effect size'] = ['NA'] * len(targets)
			if tw == 'False':
				posthoc_dfs = posthoc_dfs.drop(['Contrast', 'T'], axis=1)
			else:
				posthoc_dfs = posthoc_dfs.drop(['T'], axis=1)
			cols = ['Target Name', 'DF', 'MS', 'SS', 'p-value', 'p-value corrected', 'measures', 'distribution',
					'test', 'statistic', 'effect size']
			stats_dfs = stats_dfs.reindex(columns=cols)

			posthoc_dfs = posthoc_dfs.drop(['Contrast'], axis=1)
			posthoc_dfs = posthoc_dfs.rename(columns={'hedges': 'effect size', 'p-corr': 'p-value corrected', 'p-unc': 'p-value',
								'p-adjust': 'correction method', 'BF10': 'Bayes factor'})
			if tw == 'False':
				cols2 = ['Target Name', 'A', 'B', 'DF', 'p-value corrected', 'p-value', 'correction method', 'Paired',
					 'Parametric', 'Test', 'effect size', 'Bayes factor']
			else:
				cols2 = ['Target Name', 'Contrast', 'Group1', 'A', 'B', 'DF', 'p-value corrected', 'p-value',
						 'correction method', 'Paired',
						 'Parametric', 'Test', 'effect size', 'Bayes factor']
			posthoc_dfs = posthoc_dfs.reindex(columns=cols2)

	return stats_dfs, posthoc_dfs


# Extract groups from sample name
def add_groups(df, tw, groups1, groups2=None, colname1=None, colname2=None):
	if tw == 'False':
		df['Group'] = df['Sample Name'].str.extract(re.compile('(' + '|'.join(groups1) + ')', re.IGNORECASE),
													expand=False).fillna('')
		df['Sample Name'] = df['Sample Name'].str.replace(re.compile('(' + '|'.join(groups1) + ')', re.IGNORECASE), '')
	else:
		df['Group1'] = df[colname1].str.extract(re.compile('(' + '|'.join(groups1) + ')', re.IGNORECASE),
													expand=False).fillna('')
		df[colname1] = df[colname1].str.replace(re.compile('(' + '|'.join(groups1) + ')', re.IGNORECASE), '')
		df['Group2'] = df[colname1].str.extract(re.compile('(' + '|'.join(groups2) + ')', re.IGNORECASE),
												expand=False).fillna('')
		df[colname2] = df[colname2].str.replace(re.compile('(' + '|'.join(groups2) + ')', re.IGNORECASE), '')

	return df
