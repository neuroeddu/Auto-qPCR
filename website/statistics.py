import pandas
import numpy as np
import pingouin as pg
from pingouin import pairwise_ttests, multicomp, ttest



def stats(quantity, data, targets, rm, posthoc):

    # prepare data from intermediate dataframe
    data = data[data['Outliers'].eq(False)]
    data = data.drop(['NormQuant'], axis=1)
    data = data.drop_duplicates(keep='first')

    # write to csv
    if quantity == 2:
        # T-Test between 2 groups
        group = data['Group']
        group = group.drop_duplicates(keep='first')
        for item in targets:
            df = data[data['Target Name'].eq(item)]
            group1 = df[df['Group'].eq(group[0])]['NormMean']
            group2 = df[df['Group'].eq(group[1])]['NormMean']
            ttest(group1, group2, paired=bool(rm))
    elif quantity >= 3:
        # ANOVA test
        anova_dfs = pandas.DataFrame()
        posthoc_dfs = pandas.DataFrame()
        if rm == 'True':
            # repeated measure anova
            aov = pandas.DataFrame(['Repeated measure anova'])
            aov = aov.append(pg.rm_anova(dv='NormMean', within='Group', subject='Target Name', data=data))
            anova_dfs = aov
            ph = pandas.DataFrame([posthoc])
            ph = ph.append(posthocs_test(posthoc, data=data))
            posthoc_dfs = ph
        else:
            # anova
            pvals=[]
            for item in targets:
                aov = pg.anova(dv='NormMean', between='Group', data=data[data['Target Name'].eq(item)], detailed=True)
                pvals.append(aov['p-unc'][0])
                aov2 = pandas.DataFrame(['anova_'+item])
                aov = aov2.append(aov)
                if anova_dfs is None:
                    anova_dfs = aov
                else:
                    anova_dfs = anova_dfs.append(aov, ignore_index=True)

                print("")
                if posthoc == 'ttest_fdr':
                    ph = pandas.DataFrame(['T-Test+FDR_'+item])
                    ph = ph.append(posthocs_test('ttest_fdr', data=data, item=item))

                    if posthoc_dfs is None:
                        posthoc_dfs = ph
                    else:
                        posthoc_dfs = posthoc_dfs.append(ph, ignore_index=True)
                    print(ph)
                    print("")

            if posthoc == 'bonferroni':
                ph = pandas.DataFrame(['Bonferroni'])
                ph = ph.append(posthocs_test('bonferroni', pvals=pvals), ignore_index=True)
                posthoc_dfs = ph
                print(ph)

    return anova_dfs, posthoc_dfs


def posthocs_test(posthoc, data=None, item=None, pvals=None):
    # Post hocs
    if posthoc == 'ttest_fdr':
        post_hocs = pairwise_ttests(data=data[data['Target Name'].eq(item)], dv='NormMean' , between='Group' ,
                                        padjust='fdr_bh')

    elif posthoc == 'bonferroni':
        reject, pvals_corr = pg.multicomp(pvals, alpha=0.05, method='bonf')
        post_hocs = pandas.DataFrame({'reject': reject,
                                      'pvals_corr': pvals_corr})
    #if posthoc == 'dunnetts_test':
    return post_hocs

# def simple_anova(group, col, data):
#     col_list = data[col].drop_duplicates(keep='first')
#     for item in col_list:
#         filename = 'anova_' + item
#         aov = pg.anova(dv='NormMean' , between=group , data=data[data[col].eq(item)] , detailed=True)
#         print(item)
#         print(aov)
#         print("")

