from cgi import test
from operator import index
import unittest
import pandas as pd
from application import AUTOqPCR, routes, stability, statistics, relative, absolute, regex_rename
import os
import re

class IntegrationTest(unittest.TestCase):
    
    def __init__(self, methodName, test_data, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_folder) -> None:
        super(IntegrationTest, self).__init__(methodName)
        self.test_data = test_data
        self.model = model 
        self.quencher = quencher 
        self.task = task 
        self.cgenes = cgenes 
        self.cutoff = cutoff
        self.max_outlier = max_outlier
        self.preserve = preserve
        self.target_order = target_order
        self.sample_order = sample_order 
        self.csample = csample
        self.colnames = colnames
        self.clean_output = os.path.join(output_folder, "clean_data.csv")
        self.output = os.path.join(output_folder, "summary_data.csv")
        self.output_w_group = os.path.join(output_folder, "summary_data_w_groups.csv")

    def test_outputs(self):
        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(self.test_data, self.model, self.quencher,
                                                                                    self.task, self.cgenes, self.cutoff,
                                                                                    self.max_outlier, self.preserve,
                                                                                    self.target_order,
                                                                                    self.sample_order, self.csample,
                                                                                    self.colnames)

        clean_data_results = pd.read_csv(self.clean_output, index_col=0)

        # summary_data_results = pd.read_csv(self.output, index_col=0)

        # summary_data_w_group_results = pd.read_csv(self.output_w_group, index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)

        # targets_sorted = targets
        # samples_sorted = samples

        # if self.target_order != '':
        #     targets_sorted = []
        #     targets_sort_names = [sorter.strip() for sorter in self.target_order.split(',')]
        
        #     for name in targets_sort_names:
        #         for target in targets:
        #             if name in target and not (target in targets_sorted):
        #                 targets_sorted.append(target)

        # if self.sample_order != '':
        #     samples_sorted = []
        #     samples_sort_names = [sorter.strip() for sorter in self.sample_order.split(',')]

        #     for name in samples_sort_names:
        #         for sample in samples:
        #             if name in sample and not (sample in samples_sorted):
        #                 samples_sorted.append(sample)
		# # filter data to only have targets and samples that are mentionned


        # summary_data = summary_data.loc[targets_sorted, slice(None), :]
        # summary_data = summary_data.loc[slice(None), samples_sorted, :]

        # pd.testing.assert_frame_equal(summary_data, summary_data_results)

        # pd.testing.assert_frame_equal(summary_data_w_group, summary_data_w_group_results)

    

    def test_absolute_duplicate(self):

        # test_data = self.read_data()

        # directory = "test_data/Figure5"
        # files = ["B2M_n.csv", "NRXN3_n.csv"]
        # genes = "B2M,NRXN3"
        # for file in files:
        #     filedata = pd.read_csv(directory+"/"+file,skip_blank_lines=True, skipinitialspace=True,engine="python",encoding="utf-8",header=2)
        #     filedata["filename"] = file
        #     filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        #     filedata = filedata.loc[:,~filedata.columns.duplicated()]
        #     test_data = test_data.append(filedata, ignore_index=True, sort=True)



        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(self.test_data, self.model, self.quencher,
                                                                                            self.task, self.cgenes, self.cutoff,
                                                                                            self.max_outlier, self.preserve,
                                                                                            self.target_order,
                                                                                            self.sample_order, self.csample,
                                                                                            self.colnames)
        #clean_data_results = pd.read_csv("test_output/clean_data1.csv", index_col=0)
        
        clean_data_results = pd.read_csv(self.clean_output, index_col=0)
        
        for index, row in clean_data.iterrows():
            if pd.isnull(row["NormQuant"]):
                interest = clean_data_results[(clean_data_results["Sample Name"] == row["Sample Name"]) & (clean_data_results["Target Name"] == row["Target Name"]) & (clean_data_results["NormQuant"].isnull())]
                if interest.shape[0] != 0:
                    self.assertEqual(interest["Outliers"].bool(), row["Outliers"])
            else:
                interest = clean_data_results[(clean_data_results["Sample Name"] == row["Sample Name"]) & (clean_data_results["Target Name"] == row["Target Name"]) & (clean_data_results["NormQuant"] == row["NormQuant"])]
                if interest.shape[0] == 1:
                    self.assertEqual(interest["Outliers"].bool(), row["Outliers"])
        # clean_data_results.drop("Treatment", inplace=True, axis=1)

        # pd.testing.assert_frame_equal(clean_data, clean_data_results)

if __name__ == '__main__':
    unittest.main()