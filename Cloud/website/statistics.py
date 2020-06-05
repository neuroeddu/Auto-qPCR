import pandas
import pingouin as pg
from pingouin import pairwise_ttests, multicomp, ttest
from scipy.stats import mannwhitneyu, kruskal
import re


def stats(model, quantity, data, targets, rm, nd, posthoc):

	if model == 'absolute':
		data = data.drop(['NormQuant'], axis=1)
		mean = 'NormMean'
	elif model == 'relative':
		data = data.drop(['rq'], axis=1)
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
			group = data['Group']
			group = group.drop_duplicates(keep='first').values.tolist()
			for item in targets:
				df = data[data['Target Name'].eq(item)]
				group1 = df[df['Group'].eq(group[0])][mean]
				group2 = df[df['Group'].eq(group[1])][mean]
				t_test = ttest(group1, group2, paired=bool(rm))
				if stats_dfs is None:
					stats_dfs = t_test
				else:
					stats_dfs = stats_dfs.append(t_test, ignore_index=True)
		elif quantity >= 3:
			# ANOVA test
			stats_dfs = pandas.DataFrame()
			posthoc_dfs = pandas.DataFrame()
			if rm == 'True':
				# repeated measure anova
				aov = pandas.DataFrame(['Repeated measure anova'])
				aov = aov.append(pg.rm_anova(dv=mean, within='Group', subject='Target Name', data=data))
				stats_dfs = aov

			else:
				# anova
				pvals=[]
				for item in targets:
					aov = pg.anova(dv=mean, between='Group', data=data[data['Target Name'].eq(item)], detailed=True)
					pvals.append(aov['p-unc'][0])
					aov2 = pandas.DataFrame(['anova_'+item])
					aov = aov2.append(aov)
					if stats_dfs is None:
						stats_dfs = aov
					else:
						stats_dfs = stats_dfs.append(aov, ignore_index=True)
					if posthoc == 'ttest_fdr':
						ph = pandas.DataFrame(['T-Test+FDR_'+item])
						ph = ph.append(pairwise_ttests(data=data[data['Target Name'].eq(item)], dv=mean, between='Group',
													   padjust='fdr_bh'))
						if posthoc_dfs is None:
							posthoc_dfs = ph
						else:
							posthoc_dfs = posthoc_dfs.append(ph, ignore_index=True)
					# bonferroni test
					else:
						ph = pandas.DataFrame(['Bonferroni'])
						reject, pvals_corr = pg.multicomp(pvals, alpha=0.05, method='bonf')
						post_hocs = pandas.DataFrame({'reject': reject, 'pvals_corr': pvals_corr}, index=range(len(pvals)))
						ph = ph.append(post_hocs)
						posthoc_dfs = ph

	# nonparametric tests for not normally distributed data
	else:
		if quantity == 2:
			# Mann-Whitney U test
			stats_dfs = pandas.DataFrame()
			posthoc_dfs = pandas.DataFrame()
			# posthoc_dfs = pandas.DataFrame()
			group = data['Group']
			group = group.drop_duplicates(keep='first').values.tolist()
			for item in targets:
				df = data[data['Target Name'].eq(item)]
				group1 = df[df['Group'].eq(group[0])][mean]
				group2 = df[df['Group'].eq(group[1])][mean]
				mwu_test = mannwhitneyu(group1, group2)
				t = pandas.DataFrame(['MannWhitney_' + item])
				mwu_test = t.append(
					pandas.DataFrame({'statistic': mwu_test.statistic , 'pvalue': mwu_test.pvalue} , index=[1]))
				if stats_dfs is None:
					stats_dfs = mwu_test
				else:
					stats_dfs = stats_dfs.append(mwu_test , ignore_index=True)

		elif quantity >= 3:
			stats_dfs = pandas.DataFrame()
			posthoc_dfs = pandas.DataFrame()
			if rm == 'True':
				# friedman test for repeated measurements
				f = pandas.DataFrame(['Friedman test'])
				f = f.append(pg.friedman(dv=mean , within='Group' , subject='Target Name' , data=data))
				stats_dfs = f
			else:
				# Kruskal-Wallis H test
				pvals = []
				group = data['Group']
				group = group.drop_duplicates(keep='first').values.tolist()
				for item in targets:
					df = data[data['Target Name'].eq(item)]
					groups=[]
					for i in range(len(group)):
						groups.append(df[df['Group'].eq(group[i])][mean])
					kw_test = kruskal(*groups)
					pvals.append(kw_test.pvalue)
					t = pandas.DataFrame(['Kruskal-Wallis_' + item])
					kw_test = t.append(pandas.DataFrame({'statistic': kw_test.statistic, 'pvalue': kw_test.pvalue}, index=[1]))
					if stats_dfs is None:
						stats_dfs = kw_test
					else:
						stats_dfs = stats_dfs.append(kw_test, ignore_index=True)
					# bonferroni test
					if posthoc == 'bonferroni':
						ph = pandas.DataFrame(['Bonferroni'])
						reject, pvals_corr = pg.multicomp(pvals, alpha=0.05, method='bonf')
						post_hocs = pandas.DataFrame({'reject': reject, 'pvals_corr': pvals_corr}, index=range(len(pvals)))
						ph = ph.append(post_hocs)
						posthoc_dfs = ph
					else:
						# FDR_BH tests
						ph = pandas.DataFrame(['FDR_BH'])
						reject, pvals_corr = pg.multicomp(pvals , alpha=0.05 , method='fdr_bh')
						post_hocs = pandas.DataFrame({'reject': reject , 'pvals_corr': pvals_corr}, index=range(len(pvals)))
						ph = ph.append(post_hocs)
						posthoc_dfs = ph

	return stats_dfs, posthoc_dfs


# Extract groups from sample name
def add_groups(df, groups):
	df['Group'] = df['Sample Name'].str.extract(re.compile('(' + '|'.join(groups) + ')', re.IGNORECASE), expand=False).fillna('')
	df['Sample Name'] = df['Sample Name'].str.replace(re.compile('(' + '|'.join(groups) + ')', re.IGNORECASE), '')

	return df
