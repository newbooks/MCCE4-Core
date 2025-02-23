#!/usr/bin/env python

"""
MCCE4 Step 1: Read PDB file and generate a protein object
"""

import logging
import argparse
from mcce.pdbio import *

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'


def parse_arguments():
    helpmsg = "MCCE4 Step 1: Read PDB file and generate a mcce protein object"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", nargs='?', default="", help="PDB file")
    parser.add_argument("-s", metavar="cutoff", default="", help="SAS cutoff, default is 0.05")
    parser.add_argument("-f", metavar="ftpl_folder", default="", help="Load from this ftpl folder")
    parser.add_argument("-r", metavar="prm", nargs="+", default=[], help="Load additional runprm files, in order")
    parser.add_argument("--no_ter", default=False, action="store_true", help="Do not make terminal ressidues")
    parser.add_argument("--no_hoh", default=False, action="store_true", help="Do not include water molecules")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    logging.info("Step 1: Read PDB file and generate a mcce protein object")
    prm = Runprm()                  # Load default runprm file
    prm.update_by_cmd(args)         # Update by command line arguments
    prm.update_by_files(args.r)     # Update by additional runprm files
    prm.dump(comment="Step 1 uses these runprm parameters") # Dump to run.prm.record

