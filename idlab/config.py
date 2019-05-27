from configparser import ConfigParser, NoSectionError, NoOptionError
from pathlib import Path
import helpers as h

CONFIG = {}

# Default configuration
def set_default_config():
    global CONFIG
    CONFIG = {
        'PROJECT' : {
            # Root project directory - defaults to the system's current working directory
            'Directory': Path.cwd(),
            # Look for nested project directories
            'Recurse': False,
            # Clean and re-process project directorie with output files
            'ReRun': False,
            # Input file extensions
            'FileExt': "CSV",
            # Output file tag
            'DoneTag': "DONE"
        },
        'FILE': {
            # List of tags to extract from the files. Format is "file_line:tag_name1, file_line:tag_name2, ...
            'Tags': "",
            # CSV header file line number
            'HeaderRow': 0,
            # Commaseparated list of colum indexes to gather
            'Columns': "",
            # File encoding
            'Encoding': "utf-8"
        },
        'MODEL': {
            # The type of quantification
            'Type': "Absolute",
            # The names of the control genes
            'ControlGenes': "",
            # The names of the control samples (for stability)
            'ControlSample': "",
            # The order in which samples will appear in data tables and graphs
            'SampleOrder': "",
            # Threshold value when looking for outliers
            'OutlierCutoff': 0.0,
            # What percentage of the data set can be purged as outliers before giving up
            'MaxOutliers': 0.0
            #
        }
    }

# Update the configuration from a global configuration file
def set_global_config(sname):
    # Search for a global configuration file in the directory of the script
    sf = Path.cwd().joinpath(sname)
    cf = Path(sf).parent.joinpath(Path(sf).stem + ".conf")
    if not cf.is_file():
        h.log(h.LOG_VERBOSE, "No global configuration file found")
        return
    h.log(h.LOG_VERBOSE, "Global configuration file: {}".format(cf))
    set_config(CONFIG, cf)

# Update a configuration from a configuration file
def set_config(cfg, cfg_file):
    p = ConfigParser(
        strict = True,
        comment_prefixes = ('#'),
        allow_no_value = False,
        inline_comment_prefixes = None,
        delimiters = ('='),
        empty_lines_in_values = False)
    # Next override is to preserve the case of the options (the default is to convert them to lowercase)
    p.optionxform = lambda option: option
    p.read(cfg_file)
    # Update the configuration from the file, preserving the option value types
    for section in cfg:
        for option, value in cfg[section].items():
            try:
                if type(value) == bool:
                    cfg[section][option] = p.getboolean(section, option)
                else:
                    cfg[section][option] = (type(value))(p.get(section, option))
            # Ignore missing sections and options
            except (NoSectionError, NoOptionError):
                pass
