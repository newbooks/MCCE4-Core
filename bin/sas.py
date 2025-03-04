#!/usr/bin/env python

"""
Compute SAS for a pdb file
"""

import argparse
from mcce._sas import *

if __name__ == "__main__":
    helpmsg = "Calculate SAS for a pdb file"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", metavar="pdb", nargs=1, default=[], help="PDB file")
    parser.add_argument("--probe", metavar="probe", type=float, default=1.4, help="PDB file")
    args = parser.parse_args()


