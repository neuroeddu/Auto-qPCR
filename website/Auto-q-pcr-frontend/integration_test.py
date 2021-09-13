import unittest
import pandas as pd
from application import AUTOqPCR, routes, stability, statistics, relative, absolute, regex_rename
import os

class IntegrationTest(unittest.TestCase):

    def test_instability(self):
        test_data = pd.DataFrame()
        filename = "INSTABILITY_example.csv"

        filedata = pd.read_csv("test_data/"+filename,skip_blank_lines=True,skipinitialspace=True,engine='python',encoding="utf-8",header=46)

        filedata['filename'] = filename
        filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        filedata = filedata.loc[:,~filedata.columns.duplicated()]
        test_data = test_data.append(filedata, ignore_index=True, sort=True)

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "instability", "",
																								 "UNKNOWN", "CHR4", 0.3,
																								 0.5, False,
																								 "CHR1,CHR4,CHR8,CHR10,CHR12,CHR17,CHR18,CHR20,CHRX",
																								 "GM25953,GM25975,GM25974,GM25952,Normal", "Normal",
																								 "")

        
        clean_data_results = pd.read_csv("test_output/genomic_instability/clean_data.csv", index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)

    def test_relative(self):
        test_data = pd.DataFrame()
        filename = "RELATIVE_example.csv"

        filedata = pd.read_csv("test_data/"+filename,skip_blank_lines=True,skipinitialspace=True,engine='python',encoding="utf-8",header=46)

        filedata['filename'] = filename
        filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        filedata = filedata.loc[:,~filedata.columns.duplicated()]
        test_data = test_data.append(filedata, ignore_index=True, sort=True)

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "relative_dCT", "",
																								 "UNKNOWN", "ACTB,GAPDH", 0.3,
																								 0.5, False,
																								 "PAX6,CAMK2A,GRIN1",
																								 "AiW002-2-D0, AiW002-2-D7,KYOU-D0,KYOU-D7", "",
																								 "")

        
        clean_data_results = pd.read_csv("test_output/Relative/deltaCT/clean_data.csv", index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)

    def test_absolute(self):
        test_data = pd.DataFrame()
        directory = "test_data/Absolute"
        for filename in os.listdir(directory):
            if filename.endswith(".csv"):
                filedata = pd.read_csv(directory+"/"+filename,skip_blank_lines=True,skipinitialspace=True,engine='python',encoding="utf-8",header=46)

                filedata['filename'] = filename
                filedata.rename(columns=regex_rename.rx_rename, inplace=True)
                filedata = filedata.loc[:,~filedata.columns.duplicated()]
                test_data = test_data.append(filedata, ignore_index=True, sort=True)

        if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
            test_data.replace('Undetermined', 40,  inplace=True)

        clean_data, summary_data, summary_data_w_group, targets, samples = AUTOqPCR.process_data(test_data, "absolute", "",
																								 "UNKNOWN", "ACTB,GAPDH", 0.3,
																								 0.5, False,
																								 "",
																								 "NCRM1-IPSC,522-266-2-IPSC,AiW001-2-IPSC,AiW002-2-IPSC,AJC001-5-IPSC,AJG001C4-IPSC,NCRM1-NPC,522-266-2-NPC,AiW001-2-NPC,AiW002-2-NPC,AJC001-5-NPC, AJG001C4-NPC,NCRM1-DA4W,522-266-2-DA4W,AiW001-2-DA4W,AiW002-2-DA4W,AJG001C4-DA4W,AJC001-5-DA4W,NCRM1-DA6W,522-266-2-DA6W,AiW001-2-DA6W,AiW002-2-DA6W,AJG001C4-DA6W,AJC001-5-DA6W", "",
																								 "")

        
        clean_data_results = pd.read_csv("test_output/Absolute/clean_data.csv", index_col=0)

        pd.testing.assert_frame_equal(clean_data, clean_data_results)



if __name__ == '__main__':
    unittest.main()

