#!/usr/bin/env python

"""
Compute SAS for a pdb file
"""

import os
import sys
import argparse
import logging
from mcce._sas import sas_pdb
from mcce.pdbio import Pdb

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

if __name__ == "__main__":
    helpmsg = "Calculate SAS for a pdb file"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", metavar="pdb", nargs=1, default=[], help="PDB file")
    parser.add_argument("--probe", metavar="probe", type=float, default=1.4, help="PDB file")
    parser.add_argument("--no_hoh", default=False, action="store_true", help="remove HOH in PDB file when calculate SAS")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    pdb_file = args.pdb_file[0]

    if os.path.exists(pdb_file):
        pdb = Pdb(pdb_file)
    else:
        print(f"File {pdb_file} does not exist.")
        sys.exit(1)
    
    # check if the pdb file is actually loaded successfully
    if not pdb.mcce_ready:
        print(f"{pdb.message}")
        sys.exit(1)

    # remove water molecules if needed
    if args.no_hoh:
        pdb.remove_hoh()

    # calculate SAS
    output_lines = sas_pdb(pdb, args.probe)
    for line in output_lines:
        print(line)


