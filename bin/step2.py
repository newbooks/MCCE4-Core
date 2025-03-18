#!/usr/bin/env python

"""
MCCE4 Step 2: Make Conformers
"""

import logging
import argparse
from mcce.pdbio import *
from mcce.main import *
from mcce.constants import *

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'


def parse_arguments():
    helpmsg = "MCCE4 Step 2: Make conformers for protein residues"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("-f", metavar="ftpl_folder", default="", help="Load from this ftpl folder")
    parser.add_argument("-r", metavar="prm", nargs="+", default=[], help="Load additional runprm files, in order")
    parser.add_argument("-l", "--level", metavar="run_level", default=1, type=int, help="Run level: 1, 2, or 3, default is 1.")
    parser.add_argument("--use_head1", default=False, action="store_true", help=f"Use {STEP1_HEAD} to overwrite run.prm and run_level to make rotamers")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)


    logging.info("Step 2: Read PDB file and generate a mcce protein object")
    
    # Get prm
    prm = Runprm()                  # Create a Runprm instance and load the default runprm file
    prm.update_by_cmd(args)         # Update parameters using command line arguments
    prm.update_by_files(args.r)     # Update parameters using additional runprm files
    prm.dump(comment="Step 2 uses these runprm parameters") # Save the parameters to run.prm.record

    # Get tpl
    tpl = Tpl()                     # Create a Tpl instance
    tpl.load_ftpl_folder(prm._FTPL_FOLDER.value) # Load the ftpl folder specified in runprm
    if os.path.isdir(USER_PARAM):
        tpl.load_ftpl_folder(USER_PARAM)  # Load user defined ftpl files
    tpl.dump(comment="Step 2 uses these tpl parameters") # Save the parameters to tpl.dat.record

    # Get protein from mccepdb step1_out.pdb
    mcce = MCCE(prm=prm, tpl=tpl)
    mcce.load_mccepdb()    # load mccepdb within mcce object as tpl is required to load the mccepdb
    logging.info(f"   Protein loaded from {STEP1_OUT}")
    rot_stat = RotStat(mcce.protein)

    # count conformers at the start
    rot_stat.count_stat(mcce.protein, step="start")
    rot_stat.write_stat(mcce.protein)

    # place missing heavy atoms
    mcce.assign_qr()
    mcce.make_connect12()
    logging.info("   Place missing heavy atoms ...")
    while True:
        if mcce.place_missing_heavy_atoms() == 0:   # place_missing_heavy_atoms() returns the number of atoms placed
            break

    # make unified rotate rules
    mcce.prepare_rotate_rules()
    # mcce.print_rotate_rules()

    mcce.protein.dump(STEP2_OUT)  # Save the protein to a pdb file