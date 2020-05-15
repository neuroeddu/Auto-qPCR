import pandas
import pingouin as pg
from pingouin import pairwise_ttests, multicomp, ttest


def stats(model, quantity, data, targets, rm, posthoc):
	if model == 'absolute':
		data = data.drop(['NormQuant'], axis=1)
		mean = 'NormMean'
	elif model == 'relative':
		mean = 'rq'

	# prepare data from intermediate dataframe
	data = data[data['Outliers'].eq(False)]
	data = data.drop_duplicates(keep='first')

	# write to csv
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
			# ph = pandas.DataFrame([posthoc])
			# ph = ph.append(posthocs_test(posthoc, data=data))
			# posthoc_dfs = ph
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
				print(stats_dfs)
				print("")
				if posthoc == 'ttest_fdr':
					ph = pandas.DataFrame(['T-Test+FDR_'+item])
					ph = ph.append(pairwise_ttests(data=data[data['Target Name'].eq(item)], dv=mean, between='Group' ,
												   padjust='fdr_bh'))

					if posthoc_dfs is None:
						posthoc_dfs = ph
					else:
						posthoc_dfs = posthoc_dfs.append(ph, ignore_index=True)
					print(ph)
					print("")

			if posthoc == 'bonferroni':
				ph = pandas.DataFrame(['Bonferroni'])
				reject , pvals_corr = pg.multicomp(pvals , alpha=0.05 , method='bonf')
				post_hocs = pandas.DataFrame({'reject': reject ,
											  'pvals_corr': pvals_corr})
				ph = ph.append(post_hocs)
				posthoc_dfs = ph
				print(ph)

	return stats_dfs, posthoc_dfs


# Extract groups from sample name
def add_groups(df, groups):
	for i_row , row in df.iterrows():
		for group in groups:
			if group in df.at[i_row , 'Sample Name']:
				df.at[i_row , 'Group'] = group
				df.at[i_row , 'Sample Name'] = df.at[i_row , 'Sample Name'].replace(group , '')

	return df