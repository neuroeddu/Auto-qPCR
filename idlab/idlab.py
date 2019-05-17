# ID Lab qPCR Analysis
# Zero Clause BSD Copyright (c) 2019 by Iveta Demirova

VERSION = "0.1.5"
QUALITY = ""

from sys import argv
from pathlib import Path
from argparse import ArgumentParser
from tkinter import Tk
from tkinter.filedialog import askdirectory

import helpers as h
import config as c
import project as p


def main():
    h.verbosity = h.LOG_DEBUG
    c.set_default_config()
    if __debug__:
        if QUALITY == "DEV":
            h.verbosity = h.LOG_DEBUG
            c.CONFIG['PROJECT']['Directory'] = Path.home().joinpath \
            ("Documents/McGill/_Research/Durcan_Lab/qPCR_program/qPCR_Program-0.1.5/idlab/Gilles_plates")
    c.set_global_config(argv[0])

    #get_cmd_args()

    root = Tk()
    root.withdraw()
    c.CONFIG['PROJECT']['Directory'] = askdirectory(title = "Please select the input directory (Containing the input data and project .config file)")
    root.destroy()

    p.process_root_directory()
    h.log(h.LOG_INFO, "\nAll done")

# Update the configuration from the command line arguments
def get_cmd_args():
    p = ArgumentParser(description = "ID Lab qPCR Analysis v" + VERSION + " " + QUALITY)
    p.add_argument(
        '-d', '--directory', nargs='?', help = "project directory", dest='Directory')
    p.add_argument(
        "-r", '--recurse', dest = "Recurse", action='store_true', help = "recurse the project directory")
    p.add_argument(
        "-v", '--verbosity', dest = "Verbosity", type = int, choices = range(0,10),
        help="verbosity [{}:Quiet, {}+:Info, {}+:Verbose, {}:Debug]".format(
            h.LOG_QUIET, h.LOG_INFO, h.LOG_VERBOSE, h.LOG_DEBUG))
    args=vars(p.parse_args())

    if args["Directory"] is not None:
        c.CONFIG['PROJECT']['Directory'] = args["Directory"]

    if args["Recurse"] and not c.CONFIG['PROJECT']['Recurse']:
        c.CONFIG['PROJECT']['Recurse'] = True

    if args["Verbosity"] is not None:
        h.verbosity = args["Verbosity"]

if __name__ == "__main__":
    main()
