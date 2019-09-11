# ID Lab qPCR Analysis
# Zero Clause BSD Copyright (c) 2019 by Iveta Demirova

VERSION = "0.1.5"
QUALITY = ""

from sys import argv
from argparse import ArgumentParser, FileType
import pandas
import numpy as np
import pingouin as pg
from pingouin import pairwise_ttests, multicomp, ttest


def main():
    # Gets arguments
    args = get_cmd_args()

    # Creates empty data frame
    data = pandas.DataFrame()

    # Reads and appends all data from input files to the data dataframe
    for item in args['file']:
        filedata = pandas.read_csv(item ,
                                   skip_blank_lines=True ,
                                   skipinitialspace=True ,
                                   engine='python' ,
                                   encoding=args['encoding'] ,
                                   header=args['header'] - 1)
        data = data.append(filedata , ignore_index=True , sort=True)


    # Runs the main process model with input arguments
    #
    data, data_summary, targets, samples = process_data(data , args['mod'] , args['cgenes'] , args['ocutoff'] ,
                                                        args['omax'] , args['csample'])
    # data.to_csv('stability_result.csv')
    stats(4, data, targets, False, 'bonferroni')
    #results.to_excel(args['output'] + "output.xlsx" , encoding=args['encoding'])
    #results.to_csv(args['output'] + "output.csv")


def process_data(data , model , cgenes , cutoff , max_outliers , csample=None):
    """This filters the data and processes the selected model, returning a list of output dataframes"""

    # Transforms certain columns from string to numeric
    cols = ['CT' , 'Quantity']
    data[cols] = data[cols].apply(pandas.to_numeric , errors='coerce')

    # Marks the Control Genes in a new column in the dataframe
    data['Control'] = data['Target Name'].apply(lambda x: True if str(x) in cgenes else False)

    # Create column 'Ignore' in dataframe to mark rows with NaN values in certain columns
    data['Ignore'] = False
    data['Outliers'] = False
    cols = ['Sample Name' , 'Target Name' , 'Task' , 'Reporter' , 'CT']
    for col in cols:
        data.loc[data[col].isnull() , 'Ignore'] = True

    # Calls the different processing models depending on the model argument
    if model == 'absolute':
        data = cleanup_outliers(data , "Quantity" , cutoff , max_outliers)
        data, data_summary, targets, samples = process_absolute(data)

    elif model == 'relative':
        data = cleanup_outliers(data , "CT" , cutoff , max_outliers)
        data, data_summary, targets, samples  = process_relative(data)

    elif model == 'stability':
        data = cleanup_outliers(data , "CT" , cutoff , max_outliers)
        data, data_summary, targets, samples = process_stability(data , csample)

    return data, data_summary, targets, samples


def process_absolute(data):

    outlier_data = data[data['Outliers'].eq(True)]
    data = data[data['Outliers'].eq(False)]
    # Calculate Mean (Endogenous Control Mean) and SSD for all Controls

    data['NormQuant'] = 0

    control_filter = (data['Control'].eq(True))

    # print("Endogenous Control Quantity Means and SSD")
    # print(data_controls_grouped)

    data_controls_quantity = data[control_filter].groupby(['Sample Name']).agg({'Quantity': 'mean'})

    # print("Combined Endogenous Control Quantity Means and SSD")
    # print(data_controls_quantity)

    # Create Normalized Quantity column
    for i , row in enumerate(data_controls_quantity.itertuples(name=None) , 1):
        name_filter = (data['Sample Name'] == row[0]) #control filter used to be here
        for j in data[name_filter].index:
            data.loc[j , 'NormQuant'] = data.loc[j , 'Quantity'] / row[1]

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
    cols = ['Sample Name', 'Target Name', 'NormQuant', 'NormMean', 'NormSD', 'NormSEM', 'Outliers']
    df = pandas.DataFrame(columns=cols)
    for item in cols:
        df[item] = data[item]

    for i_row, row in df.iterrows():
        if 'IPSC' in df.at[i_row, 'Sample Name']:
            df.at[i_row, 'Group'] = 'IPSC'
            df.at[i_row, 'Sample Name'] = df.at[i_row, 'Sample Name'].replace('-IPSC', '')
        elif 'NPC' in df.at[i_row, 'Sample Name']:
            df.at[i_row, 'Group'] = 'NPC'
            df.at[i_row, 'Sample Name'] = df.at[i_row , 'Sample Name'].replace('-NPC', '')
        elif 'DA4W' in df.at[i_row, 'Sample Name']:
            df.at[i_row, 'Group'] = 'DA4W'
            df.at[i_row, 'Sample Name'] = df.at[i_row , 'Sample Name'].replace('-DA4W' , '')
        elif 'DA6W' in df.at[i_row, 'Sample Name']:
            df.at[i_row, 'Group'] = 'DA6W'
            df.at[i_row, 'Sample Name'] = df.at[i_row , 'Sample Name'].replace('-DA6W' , '')

    data_output_summary = data.groupby(['Target Name' , 'Sample Name']).agg(
        {'NormQuant': [np.size , 'mean' , 'std'] , 'NormSEM': 'mean'})

    return df, data_output_summary, targets, samples


def process_relative(data):

    outlier_data = data[data['Outliers'].eq(True)]
    data = data[data['Outliers'].eq(False)]

    # Calculate CT Mean (Endogenous Control Mean) and SSD for all Controls

    control_filter = (data['Control'].eq(True))
    data_controls_grouped = data[control_filter].groupby(['Target Name' , 'Sample Name']).agg(
        {'CT': [np.size , 'mean' , 'std']})

    # print("Endogenous Control CT Means and SSD")
    # print(data_controls_grouped)

    data_controls_ct = data[control_filter].groupby(['Sample Name']).agg({'CT': 'mean'})

    # print("Combined Endogenous Control CT Means and SSD")
    # print(data_controls_ct)

    # Create rq column
    for i , row in enumerate(data_controls_ct.itertuples(name=None) , 1):
        name_filter = (data['Sample Name'] == row[0])
        for j in data[name_filter].index:
            data.loc[j , 'rq'] = np.power(2 , -(data.loc[j , 'CT'] - row[1]))

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
            mean_sem_result[target][sample] = (mean , sdt_dev , std_err)
    for i_row , row in data.iterrows():
        if data.at[i_row , 'Sample Name'] in samples and data.at[i_row , 'Sample Name'] in samples and data.at[
            i_row , 'Target Name'] in mean_sem_result and data.at[i_row , 'Sample Name'] in mean_sem_result[
            data.at[i_row , 'Target Name']]:
            data.at[i_row , 'rq'] = mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][0]
            data.at[i_row , 'rqSD'] = mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][1]
            data.at[i_row , 'rqSEM'] = mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][
                2]
    # Making the intermediate dataframe
    data = data.append(outlier_data)
    cols = ['Sample Name' , 'Target Name' , 'rq' , 'rqSD' , 'rqSEM' , 'Outliers']
    df = pandas.DataFrame(columns=cols)
    for item in cols:
        df[item] = data[item]

    data_output_summary = data.groupby(['Target Name' , 'Sample Name']).agg(
        {'rq': [np.size , 'mean'] , 'rqSD': 'mean' , 'rqSEM': 'mean'})

    return df, data_output_summary, targets, samples


def process_stability(data , csample):
    outlier_data = data[data['Outliers'].eq(True)]
    data = data[data['Outliers'].eq(False)]

    # Calculate CT Mean (Endogenous Control Mean) and SSD for all Controls

    control_filter = (data['Control'].eq(True))
    data_controls_grouped = data[control_filter].groupby(['Target Name' , 'Sample Name']).agg(
        {'CT': [np.size , 'mean' , 'std']})

    # print("Endogenous Control CT Means and SSD")
    # print(data_controls_grouped)

    data_controls_ct = data[control_filter].groupby(['Sample Name']).agg({'CT': 'mean'})

    # print("Combined Endogenous Control CT Means and SSD")
    # print(data_controls_ct)

    # Create deltaCT column
    for i , row in enumerate(data_controls_ct.itertuples(name=None) , 1):
        name_filter =  (data['Sample Name'] == row[0])
        for j in data[name_filter].index:
            data.loc[j , 'deltaCT'] = data.loc[j , 'CT'] - row[1]

    # Mark the Control Samples
    data['ControlSample'] = data['Sample Name'].apply(lambda x: True if x in csample else False)
    filter_sample = (data['ControlSample'].eq(True))
    data_sampled = data[filter_sample].groupby(['Target Name']).agg({('deltaCT'): 'mean'})

    # print("Mean Control Sample Delta CT")
    # print(data_sampled)

    for i , row in enumerate(data_sampled.itertuples(name=None) , 1):
        target_filter = (data['Target Name'] == row[0])
        for j in data[target_filter].index:
            data.loc[j , 'rq'] = np.power(2 , -(data.loc[j , 'deltaCT'] - row[1]))

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
            mean_sem_result[target][sample] = (mean , sdt_dev , std_err)
    for i_row , row in data.iterrows():
        if data.at[i_row , 'Sample Name'] in samples and data.at[i_row , 'Sample Name'] in samples and data.at[
            i_row , 'Target Name'] in mean_sem_result and data.at[i_row , 'Sample Name'] in mean_sem_result[
            data.at[i_row , 'Target Name']]:
            data.at[i_row , 'rq'] = mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][0]
            data.at[i_row , 'rqSD'] = mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][1]
            data.at[i_row , 'rqSEM'] = mean_sem_result[data.at[i_row , 'Target Name']][data.at[i_row , 'Sample Name']][
                2]

    # Making the intermediate dataframe
    data = data.append(outlier_data)
    cols = ['Sample Name' , 'Target Name' , 'rq' , 'rqSD' , 'rqSEM' , 'Outliers']
    df = pandas.DataFrame(columns=cols)
    for item in cols:
        df[item] = data[item]

    data_output_summary = data.groupby(['Target Name' , 'Sample Name']).agg(
        {'rq': [np.size , 'mean'] , 'rqSD': 'mean' , 'rqSEM': 'mean'})

    return df, data_output_summary, targets, samples


def cleanup_outliers(d , feature , cutoff , max_outliers):
    """Function to remove outliers based on cutoff and maximum number of outliers,
    by removing the furthest data point in each group when the standard deviation
    is higher than the cutoff"""

    # Calculate SSD for all sample groups
    f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN')
    d1 = d[f].groupby(['Sample Name' , 'Target Name']).agg({'CT': ['std']})
    f = (d1['CT']['std'] > cutoff)
    d2 = d1[f]
    if not d2.empty:
        # Mark all outliers
        for i , row in enumerate(d2.itertuples(name=None) , 1):
            f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
                & (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
            dx_idx = d[f].index
            group_size = len(dx_idx)
            min_size = round(group_size * (1 - max_outliers))
            size = group_size
            if min_size < 2:
                min_size = 2
                h.log(5 , 'Minimum group size must be equal or greater than 2')
            while True:
                f = (d['Ignore'].eq(False)) & (d['Task'] == 'UNKNOWN') \
                    & (d['Sample Name'] == row[0][0]) & (d['Target Name'] == row[0][1])
                dx = d[f].copy()
                dxg = d[f].groupby(['Sample Name' , 'Target Name']).agg({feature: [np.size , 'std' , 'mean']})
                if dxg[feature]['std'].iloc[0] <= cutoff:
                    # CT std is under the threshold
                    break
                # Will ignore one or all measurements
                size -= 1
                if size < min_size:
                    # Ignore the entire group of measurements
                    # for j in dx_idx:
                    #    d['Ignore'].loc[j] = True
                    break
                # Will remove the measurement which is furthest from the mean
                dx['Distance'] = (dx[feature] - dxg[feature]['mean'].iloc[0]) ** 2
                j = dx.sort_values(by='Distance' , ascending=False).index[0]
                d['Outliers'].loc[j] = True
                d['Ignore'].loc[j] = True

    return (d[(d['Ignore'].eq(False))])


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



def get_cmd_args():
    """This sets up the ArgumentParser and returns the list of arguments"""

    # Creates the Argument Parser
    parser = ArgumentParser(description="ID Lab qPCR Analysis v" + VERSION + " " + QUALITY)

    # Adds the input file argument
    parser.add_argument('-f' , '--file' ,
                        nargs='+' ,
                        type=FileType('r') ,
                        required=True)

    # Adds the output directory
    parser.add_argument('-o' , '--output' ,
                        required=True)

    # Adds the model argument, to select between the three models
    parser.add_argument('-m' , '--mod' , '--model' ,
                        nargs='?' ,
                        choices=['relative' , 'absolute' , 'stability'] ,
                        required=True)

    # Adds the control genes argument, taking a list of gene names
    parser.add_argument('-cg' , '--cgenes' , '--controlgenes' ,
                        nargs='+' ,
                        required=True)

    # Adds the optional control sample argument for the stability model, taking a list of sample names
    parser.add_argument('-cs' , '--csample' , '--controlsamples' ,
                        nargs='*')

    # Adds optional outlier cutoff
    parser.add_argument('-oc' , '--ocutoff' ,
                        type=float ,
                        default=0.3)

    # Adds optional max outliers
    parser.add_argument('-om' , '--omax' ,
                        type=float ,
                        default=0.5)

    # Adds optional encoding
    parser.add_argument('-e' , '--encoding' ,
                        default='ISO-8859-1')

    # Adds optional header size
    parser.add_argument('-hd' , '--header' ,
                        default=47)

    return vars(parser.parse_args("-f /Users/admin/Documents/Gilles_data/Gilles_data_combined.csv -o /Users/admin/Documents/Gilles_data/ -m absolute -cg GAPDH ACTB".split()))

#-f /Users/admin/Documents/GitHub/Auto-qPCR/Auto-qPCR-program/example_stability/2019-05-23_133411-control_and_ko-line.csv -o /Users/admin/Documents/GitHub/Auto-qPCR/Auto-qPCR-program/example_stability/ -m stability -cg ACTB GAPDH -cs H9".split()))
#-f /Users/admin/Documents/Gilles_data/Gilles_data_combined.csv -o /Users/admin/Documents/Gilles_data/ -m absolute -cg GAPDH ACTB".split())

#-f /Users/admin/Documents/GitHub/Auto-qPCR/Auto -qPCR-program/example_relative/data.csv -o /Users/admin/Documents/GitHub/Auto-qPCR/Auto-qPCR-program/example_relative/ -m relative -cg NANOG C-MYC ZFP42

#return vars(parser.parse_args("-f /Users/admin/Documents/GitHub/Auto-qPCR/Auto-qPCR-program/example_relative/data.csv -o /Users/admin/Documents/GitHub/Auto-qPCR/Auto-qPCR-program/example_relative/ -m relative -cg NANOG C-MYC ZFP42".split()))



# -f C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_stability/2019-05-23_133411-controlandko-line.csv -o C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_relative/ -m stability -cg GAPDH ACTB -cs H9
# -f C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_relative/data.csv -o C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_relative/ -m relative -cg GAPDH1 GAPDH2 ACTB
# -f C:/Users/eddie/Documents/idlab/test/Iva/2018-09-14GAPDHP2X3.csv C:/Users/eddie/Documents/idlab/test/Iva/2018-09-18OPRD1bACTIN.csv C:/Users/eddie/Documents/idlab/test/Iva/2018-09-20PIEZO2.csv C:/Users/eddie/Documents/idlab/test/Iva/2018-10-02RET.csv -o C:/Users/eddie/Documents/Auto-qPCR/Auto-qPCR-program/example_relative/ -m absolute -cg GAPDH ACTB BACT
# -f /Users/admin/Documents/GitHub/Auto-qPCR/Auto-qPCR-program/example_relative/data.csv -o /Users/admin/Documents/GitHub/Auto-qPCR/Auto-qPCR-program/example_relative/ -m relative -cg NANOG C-MYC ZFP42



if __name__ == "__main__":
    main()
