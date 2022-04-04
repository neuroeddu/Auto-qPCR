from cgi import test
from operator import index
import unittest
import pandas as pd
from application import AUTOqPCR, routes, stability, statistics, relative, absolute, regex_rename
import os
import re

class IntegrationTest(unittest.TestCase):
    

    def __init__(self, methodName, input_file, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_file) -> None:
        super(IntegrationTest, self).__init__(methodName)
        self.input_file = input_file
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
        self.output_file = output_file

    def test_instability(self):
        test_data = pd.DataFrame()

        #filename = "INSTABILITY_example.csv"

        #filedata = pd.read_csv("test_data/"+filename,skip_blank_lines=True,skipinitialspace=True,engine='python',encoding="utf-8",header=46)

        #filedata['filename'] = filename
        
        for file in self.input_file:
            filedata = pd.read_csv(file, skip_blank_lines=True, skipinitialspace=True, engine='python', encoding='utf-8', header=46)   
            filedata['filename'] = file
            filedata.rename(columns=regex_rename.rx_rename, inplace=True)
            filedata = filedata.loc[:,~filedata.columns.duplicated()]
            test_data = test_data.append(filedata, ignore_index=True, sort=True)

        #filedata = pd.read_csv(self.input_file, skip_blank_lines=True, skipinitialspace=True, engine='python', encoding='utf-8', header=46)

        # filedata['filename'] = self.input_file

        # filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        # filedata = filedata.loc[:,~filedata.columns.duplicated()]
        # test_data = test_data.append(filedata, ignore_index=True, sort=True)

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        # clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "instability", "",
		# 																						 "UNKNOWN", "CHR4", 0.3,
		# 																						 0.5, False,
		# 																						 "CHR1,CHR4,CHR8,CHR10,CHR12,CHR17,CHR18,CHR20,CHRX",
		# 																						 "GM25953,GM25975,GM25974,GM25952,Normal", "Normal",
		# 																						 "")

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, self.model, self.quencher,
                                                                                            self.task, self.cgenes, self.cutoff,
                                                                                            self.max_outlier, self.preserve,
                                                                                            self.target_order,
                                                                                            self.sample_order, self.csample,
                                                                                            self.colnames)

        clean_data_results = pd.read_csv(self.output_file, index_col=0)

        # clean_data_results = pd.read_csv("test_output/genomic_instability/clean_data.csv", index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)

    def test_relative(self):
        test_data = pd.DataFrame()

        for file in self.input_file:
            filedata = pd.read_csv(file, skip_blank_lines=True, skipinitialspace=True, engine='python', encoding='utf-8', header=46)   
            filedata['filename'] = file
            filedata.rename(columns=regex_rename.rx_rename, inplace=True)
            filedata = filedata.loc[:,~filedata.columns.duplicated()]
            test_data = test_data.append(filedata, ignore_index=True, sort=True)
        # filename = "RELATIVE_example.csv"

        # filedata = pd.read_csv("test_data/"+filename,skip_blank_lines=True,skipinitialspace=True,engine='python',encoding="utf-8",header=46)

        # filedata['filename'] = filename
        # filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        # filedata = filedata.loc[:,~filedata.columns.duplicated()]
        # test_data = test_data.append(filedata, ignore_index=True, sort=True)

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, self.model, self.quencher,
                                                                                            self.task, self.cgenes, self.cutoff,
                                                                                            self.max_outlier, self.preserve,
                                                                                            self.target_order,
                                                                                            self.sample_order, self.csample,
                                                                                            self.colnames)

        # clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "relative_dCT", "",
		# 																						 "UNKNOWN", "ACTB,GAPDH", 0.3,
		# 																						 0.5, False,
		# 																						 "PAX6,CAMK2A,GRIN1",
		# 																						 "AiW002-2-D0, AiW002-2-D7,KYOU-D0,KYOU-D7", "",
		# 																						 "")

        
        #clean_data_results = pd.read_csv("test_output/Relative/deltaCT/clean_data.csv", index_col=0)
        clean_data_results = pd.read_csv(self.output_file, index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)
        #pd.testing.assert_frame_equal(clean_data, clean_data_results)

    def test_absolute(self):
        test_data = pd.DataFrame()

        for file in self.input_file:
            filedata = pd.read_csv(file, skip_blank_lines=True, skipinitialspace=True, engine='python', encoding='utf-8', header=46)   
            filedata['filename'] = file
            filedata.rename(columns=regex_rename.rx_rename, inplace=True)
            filedata = filedata.loc[:,~filedata.columns.duplicated()]
            test_data = test_data.append(filedata, ignore_index=True, sort=True)
        # directory = "test_data/Absolute"
        # for filename in os.listdir(directory):
        #     if filename.endswith(".csv"):
        #         filedata = pd.read_csv(directory+"/"+filename,skip_blank_lines=True,skipinitialspace=True,engine='python',encoding="utf-8",header=46)

        #         filedata['filename'] = filename
        #         filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        #         filedata = filedata.loc[:,~filedata.columns.duplicated()]
        #         test_data = test_data.append(filedata, ignore_index=True, sort=True)

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, self.model, self.quencher,
                                                                                            self.task, self.cgenes, self.cutoff,
                                                                                            self.max_outlier, self.preserve,
                                                                                            self.target_order,
                                                                                            self.sample_order, self.csample,
                                                                                            self.colnames)

        # clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "absolute", "",
		# 																						 "UNKNOWN", "ACTB,GAPDH", 0.3,
		# 																						 0.5, False,
		# 																						 "",
		# 																						 "NCRM1-IPSC,522-266-2-IPSC,AiW001-2-IPSC,AiW002-2-IPSC,AJC001-5-IPSC,AJG001C4-IPSC,NCRM1-NPC,522-266-2-NPC,AiW001-2-NPC,AiW002-2-NPC,AJC001-5-NPC, AJG001C4-NPC,NCRM1-DA4W,522-266-2-DA4W,AiW001-2-DA4W,AiW002-2-DA4W,AJG001C4-DA4W,AJC001-5-DA4W,NCRM1-DA6W,522-266-2-DA6W,AiW001-2-DA6W,AiW002-2-DA6W,AJG001C4-DA6W,AJC001-5-DA6W", "",
		# 																						 "")

        
        clean_data_results = pd.read_csv(self.output_file, index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)

        #clean_data_results = pd.read_csv("test_output/Absolute/clean_data.csv", index_col=0)

        #pd.testing.assert_frame_equal(clean_data, clean_data_results)

    def test_delta_relative(self):
        test_data = pd.DataFrame()
        #filename = "RELATIVE_example.csv"

        #filedata = pd.read_csv("test_data/"+filename,skip_blank_lines=True,skipinitialspace=True,engine='python',encoding="utf-8",header=46)

        #filedata['filename'] = filename
        
        for file in self.input_file:
            filedata = pd.read_csv(file, skip_blank_lines=True, skipinitialspace=True, engine='python', encoding='utf-8', header=46)
            filedata['filename'] = file
            filedata.rename(columns=regex_rename.rx_rename, inplace=True)
            filedata = filedata.loc[:,~filedata.columns.duplicated()]
            test_data = test_data.append(filedata, ignore_index=True, sort=True)

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        #clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "relative_ddCT", "", "UNKNOWN", "ACTB,GAPDH", 0.3, 0.5, False, "PAX6,CAMK2A,GRIN1", "AiW002-2-D0, AiW002-2-D7,KYOU-D0,KYOU-D7", "AiW002-2-D0", "")

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, self.model, self.quencher,
                                                                                            self.task, self.cgenes, self.cutoff,
                                                                                            self.max_outlier, self.preserve,
                                                                                            self.target_order,
                                                                                            self.sample_order, self.csample,
                                                                                            self.colnames)

        #clean_data_results = pd.read_csv("test_output/Relative/deltadelatCT/clean_data.csv", index_col=0)

        clean_data_results = pd.read_csv(self.output_file, index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)

    

    def test_absolute_duplicate(self):
        test_data = pd.DataFrame()

        for file in self.input_file:
            filedata = pd.read_csv(file, skip_blank_lines=True, skipinitialspace=True, engine='python', encoding='utf-8', header=46)
            filedata['filename'] = file
            filedata.rename(columns=regex_rename.rx_rename, inplace=True)
            filedata = filedata.loc[:,~filedata.columns.duplicated()]
            test_data = test_data.append(filedata, ignore_index=True, sort=True)
        # directory = "test_data/Figure5"
        # files = ["B2M_n.csv", "NRXN3_n.csv"]
        # genes = "B2M,NRXN3"
        # for file in files:
        #     filedata = pd.read_csv(directory+"/"+file,skip_blank_lines=True, skipinitialspace=True,engine="python",encoding="utf-8",header=2)
        #     filedata["filename"] = file
        #     filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        #     filedata = filedata.loc[:,~filedata.columns.duplicated()]
        #     test_data = test_data.append(filedata, ignore_index=True, sort=True)

        genes = [genes.strip() for genes in genes.split(",")]
        test_data['Target Name'] = test_data['filename'].str.extract(re.compile('(' + '|'.join(genes) + ')', re.IGNORECASE),
																   expand=False).fillna('')

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, self.model, self.quencher,
                                                                                            self.task, self.cgenes, self.cutoff,
                                                                                            self.max_outlier, self.preserve,
                                                                                            self.target_order,
                                                                                            self.sample_order, self.csample,
                                                                                            self.colnames)

        # clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "absolute", "TMR", "sample", "B2M", 0.3, 0.5,True,"", "B4bisNST,B4bisGP,B4bisSN,B6NST,R6 NST,R6 NST,R6 GP,R6 SN,V3 NST,V3 GP,V3 SN, V4 NST,V4 GP,V4 SN,R5bis NST,R5bis GP,R5bis SN,R6bisNST,R6bisGP,R6bisSN,R8bisNST,R8bisGP, R8bisSN,V2NST,V2GP,V2SN,V8NST,V8GP,V8SN", "", "")

        #clean_data_results = pd.read_csv("test_output/clean_data1.csv", index_col=0)
        
        clean_data_results = pd.read_csv(self.output_file, index_col=0)
        
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