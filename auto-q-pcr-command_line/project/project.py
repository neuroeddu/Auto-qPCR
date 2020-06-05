import os
import re
from pathlib import Path
from copy import deepcopy
import pandas as pd
import helpers as h
import config as c
import model as m
import model_stability as ms
import model_relative as mr


INIT = {
    'project': {'local_config': None, 'data': pd.DataFrame()}
}

PROJECTS = {}

def process_root_directory():
    root = Path(c.CONFIG['PROJECT']['Directory'])
    if not root.is_dir():
        h.log(h.LOG_QUIET, 'ABORTING: The specified root directory was not found: {}'.format(root))
        quit
    # Walk the root directory
    for wdir, wdirs, wfiles in os.walk(root, topdown=True):
        stop_recurse = process_directory(wdir, wfiles)
        if stop_recurse:
            wdirs.clear()

    n = len(PROJECTS)
    h.log(h.LOG_INFO, "{} project{} found and processed".format(
        "No" if n == 0 else n, '' if n == 1 else 's'))

def process_directory(wdir, wfiles):
    wd = Path(wdir)
    h.log(h.LOG_VERBOSE, "INFO: In folder: {}".format(wd))

    # Search for a local configuration file in the current directory
    cf = wd.joinpath("project_config.conf")
    if cf.is_file():
        h.log(h.LOG_VERBOSE, "INFO: Local configuration file found")
        cfg = deepcopy(c.CONFIG)
        c.set_config(cfg, cf)
    else:
        # Use the global configuration
        cf = None
        cfg = c.CONFIG

    # Look for files with extensions specified in the configuration
    extns = set([x.strip().lower() for x in cfg['PROJECT']['FileExt'].split(',')])
    files = []
    for ext in extns:
        rx = re.compile(str.replace(h.RXP['FileExt'],"@1", ext), re.IGNORECASE)
        files += list(filter(rx.findall, wfiles))

    if len(files) > 0:
        load_files(wdir, cfg, files)

    if wdir in PROJECTS:
        PROJECTS[wdir]['local_config'] = (cf is not None)
        process_project(wdir, PROJECTS[wdir], cfg)
        
    # Should we stop walking further down this path
    return not cfg['PROJECT']['Recurse']

def load_files(wdir, cfg, files):
    # Create new project record
    PROJECTS[wdir] = deepcopy(INIT['project'])

    s = set()
    rx = re.compile(h.RXP['IntRange'], re.IGNORECASE)
    for x in cfg['FILE']['Columns'].split(','):
        x = x.strip()
        if rx.match(x):
            t = x.split('-')
            for i in range(int(t[0].rstrip()), int(t[1].lstrip())+1):
                s.add(i-1)
        else:
            s.add(int(x)-1)
    cols = list(s)

    tags = {}
    tags_idx = set()
    for x in cfg['FILE']['Tags'].split(','):
        t = x.split(':')
        a = int(t[0].strip())-1
        tags[str(a)] = t[1].strip()
        tags_idx.add(a)

    for fname in files:
        h.log(h.LOG_VERBOSE, "INFO: Loading file: {}".format(fname))
        fn = Path(wdir).joinpath(fname)

        # Load csv file data
        if len(cols) > 0:
            data = pd.read_csv(fn,
                usecols = cols,   
                skip_blank_lines = True,
                skipinitialspace = True,
                engine = 'python',
                encoding = cfg['FILE']['Encoding'],
                header = cfg['FILE']['HeaderRow'] - 1)
        else:
            data = pd.read_csv(fn,
                skip_blank_lines = True,
                skipinitialspace = True,
                engine = 'python',
                encoding = cfg['FILE']['Encoding'],
                header = cfg['FILE']['HeaderRow'] - 1)

        if not data.empty:
            # Get additional data from file
            if len(tags_idx) > 0:
                with open(fn, 'r', encoding = cfg['FILE']['Encoding']) as f:
                    for i, line in enumerate(f):
                        if i >= cfg['FILE']['HeaderRow']:
                            break
                        if i in tags_idx:
                            data.insert(0, tags[str(i)], line.split(',')[1].strip())
            # Append to the project's data
            PROJECTS[wdir]['data'] = PROJECTS[wdir]['data'].append(data, ignore_index = True)

    # Delete the project record if no data
    if  PROJECTS[wdir]['data'].empty:
        del PROJECTS[wdir]

def process_project(wdir, p, cfg):
    # if h.verbosity == h.LOG_DEBUG:
        # print(p)
    data = p['data']

    # Branches based on type of quantification specified in the .config file
    if cfg['MODEL']['Type'] == 'Stability':
        ms.run_model(wdir, data, cfg)
    if cfg['MODEL']['Type'] == 'Relative':
        mr.run_model(wdir, data, cfg)
    if cfg['MODEL']['Type'] == 'Absolute':
        m.run_model(wdir, data, cfg)
