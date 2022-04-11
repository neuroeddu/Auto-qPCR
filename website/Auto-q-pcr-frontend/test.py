import os, re
import sys, getopt
import unittest
import pandas as pd
from application import regex_rename
from integration_test import IntegrationTest

boolean_converter = {"True": True, "False": False}

def main(argv):
    input_file = ""
    output_folder = ""
    try:
        opts, arg = getopt.getopt(argv, "hi:o:", ["ifile=", "odir="])
    except getopt.GetoptError:
        print("test.py -i <input_file> -o <output_folder>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("test.py -i <input_file> -o <output_folder>")
            sys.exit()
        elif opt in ("-o", "--odir"):
            output_folder = arg
        elif opt in ("-i", "--ifile"):
            input_file = arg
    input_file = input_file.split(",")
    model = ""
    quencher = ""
    task = ""
    cgenes = ""
    genes = ""
    cutoff = ""
    max_outlier = ""
    preserve = ""
    target_order = ""
    sample_order = ""
    csample = ""
    colnames = ""
    f = open(os.path.join(output_folder, "log.txt"), "r")
    line = f.readline()
    while line:
        rline = line.strip('\n')
        if "Model" in rline and len(rline.split(":")) > 1:
            model = rline.split(":")[-1].strip()
        elif "Quencher" in rline and len(rline.split(":")) > 1:
            quencher = rline.split(":")[-1].strip()
        elif "Task" in rline and len(rline.split(":")) > 1:
            task = rline.split(":")[-1].strip()
        elif "Endogenous control genes" in rline and len(rline.split(":")) > 1:
            cgenes = rline.split(":")[-1].strip()
        elif "Cut-off" in rline and len(rline.split(":")) > 1:
            cutoff = float(rline.split(":")[-1].strip())
        elif "Maximum Outliers" in rline and len(rline.split(":")) > 1:
            max_outlier = float(rline.split(":")[-1].strip())
        elif "Preserve highly variable replicates" in rline and len(rline.split(":")) > 1:
            preserve = boolean_converter[rline.split(":")[-1].strip()]
        elif "Target Order" in rline and len(rline.split(":")) > 1:
            target_order = rline.split(":")[-1].strip()
        elif "Sample Order" in rline and len(rline.split(":")) > 1:
            sample_order = rline.split(":")[-1].strip()
        elif "Control Sample" in rline and len(rline.split(":")) > 1:
            csample = rline.split(":")[-1].strip()
        elif "Additional column names" in rline and len(rline.split(":")) > 1:
            colnames = rline.split(":")[-1].strip()
        elif "Gene names if they are included in file names" in rline and len(rline.split(":")) > 1:
            genes = rline.split(":")[-1].strip()
        line = f.readline()
    f.close()
    t = unittest.TestSuite()
    test_data = read_data(input_file, genes)
    if model == "absolute" and preserve:
        t.addTest(IntegrationTest('test_absolute_duplicate', test_data, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_folder))
    else:
        t.addTest(IntegrationTest('test_outputs', test_data, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_folder))
    unittest.TextTestRunner().run(t)

def read_data(input_file, genes):

    test_data = pd.DataFrame()

    for file in input_file:
        filedata = pd.read_csv(file, skip_blank_lines=True, skipinitialspace=True, engine='python', encoding='utf-8', header=46)   
        filedata['filename'] = os.path.basename(file)
        filedata.rename(columns=regex_rename.rx_rename, inplace=True)
        filedata = filedata.loc[:,~filedata.columns.duplicated()]
        test_data = test_data.append(filedata, ignore_index=True, sort=True)

    if genes != '':
        genes = [genes.strip() for genes in genes.split(",")]
        test_data['Target Name'] = test_data['filename'].str.extract(re.compile('(' + '|'.join(genes) + ')', re.IGNORECASE),
																   expand=False).fillna('')

    if len(test_data[test_data['CT'].astype(str).str.contains('Undetermined', na = False)]) > 0:
        test_data.replace('Undetermined', 40,  inplace=True)

    return test_data


if __name__ == "__main__":
    main(sys.argv[1:])