import sys, getopt
import unittest
from integration_test import IntegrationTest

boolean_converter = {"True": True, "False": False}

def main(argv):
    input_file = ""
    log_file = ""
    output_file = ""
    try:
        opts, arg = getopt.getopt(argv, "hi:l:o:", ["ifile=", "lfile=", "ofile="])
    except getopt.GetoptError:
        print("test.py -i <input_file> -l <log_file> -o <output_file>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("test.py -i <input_file> -l <log_file> -o <output_file>")
            sys.exit()
        elif opt in ("-o", "--ofile"):
            output_file = arg
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-l", "--lfile"):
            log_file = arg
    input_file = input_file.split(",")
    model = ""
    quencher = ""
    task = ""
    cgenes = ""
    cutoff = ""
    max_outlier = ""
    preserve = ""
    target_order = ""
    sample_order = ""
    csample = ""
    colnames = ""
    f = open(log_file, "r")
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
        line = f.readline()
    f.close()
    t = unittest.TestSuite()
    if model == "absolute":
        if preserve:
            t.addTest(IntegrationTest('test_absolute_duplicate', input_file, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_file))
        else:
            t.addTest(IntegrationTest('test_absolute', input_file, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_file))
    elif model == "instability":
        t.addTest(IntegrationTest('test_instability', input_file, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_file))
    elif model == "relative_ddCT":
        t.addTest(IntegrationTest('test_delta_relative', input_file, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_file))
    elif model == "relative_dCT":
        t.addTest(IntegrationTest('test_relative', input_file, model, quencher, task, cgenes, cutoff, max_outlier, preserve, target_order, sample_order, csample, colnames, output_file))
    unittest.TextTestRunner().run(t)


if __name__ == "__main__":
    main(sys.argv[1:])