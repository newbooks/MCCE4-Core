#!/usr/bin/env python

"""
MCCE4 Step 1: Read PDB file and generate a protein object
"""

import logging
import argparse

def parse_arguments():
    helpmsg = "MCCE4 Step 1: Read PDB file and generate a protein object"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", help="PDB file")
    parser.add_argument("-s", metavar="cutoff", type=float, default=0.05, help="SAS cutoff, default is 0.05")
    parser.add_argument("-f", metavar="ftpl_folder", default="", help="Load from this ftpl folder")
    parser.add_argument("-r", metavar="prm", nargs="+", default=[], help="Load additional runprm files, in order")
    parser.add_argument("--no_ter", default=False, action="store_true", help="Do not make terminal ressidues")
    parser.add_argument("--no_hoh", default=False, action="store_true", help="Do not include water molecules")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

