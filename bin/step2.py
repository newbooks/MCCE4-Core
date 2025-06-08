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
    parser.add_argument("--writepdb", default=False, action="store_true", help=f"Write GA microstates to mccepdb under {GA_OUTPUT_FOLDER}")
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
    if os.path.isfile(NEW_FTPL):
        tpl.load_ftpl_file(NEW_FTPL)    # Load new.ftpl
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

    # make swap conformers
    logging.info("   Make swap conformers ...")
    mcce.propogate_swap()
    rot_stat.count_stat(mcce.protein, step="swap")
    rot_stat.write_stat(mcce.protein)

    # print residue id, conformer id, and conftype for each conformer
    # for res in mcce.protein.residues:
    #     for conf in res.conformers:
    #         print(f"   {res.resid} {res.resname} {conf.confid} {conf.conftype}")


    # expand conformers by conftypes
    logging.info("   Propogate conformers with conformer types and add H atoms ...")
    # propogate conftypes
    mcce.propogate_conftypes()
    rot_stat.count_stat(mcce.protein, step="type")
    rot_stat.write_stat(mcce.protein)

    # make unified rotate rules
    logging.info("   Prepare rotamer making rules ...")
    mcce.prepare_rotate_rules()
    mcce.apply_rotate_rules()
    # mcce.print_rotate_rules()

    # explore conformers with Genetic Algorithm
    logging.info("   Explore conformers with Genetic Algorithm ...")
    logging.info(f"      This may take a while, check {GA_PROGRESS} for GA progress.")
    mcce.ga_optimize(writepdb=args.writepdb)
    rot_stat.count_stat(mcce.protein, step="ga")
    rot_stat.write_stat(mcce.protein)
    logging.info("   Done with Genetic Algorithm conformer search.")


    mcce.protein.dump(STEP2_OUT)  # Save the protein to a pdb file

    # print summary
    
    print("====================================================================")
    print("Step 2 completed. Here is run summary:")
    print("--------------------------------------------------------------------")    
    print(f"  {STEP2_OUT:<16s}: Protein object in mccepdb file, used in step 3")
    print(f"  {GA_PROGRESS:<16s}: GA progress log, information only")
    if args.writepdb:
        print(f"  {GA_OUTPUT_FOLDER + '/':<16s}: GA selected individuals in mccepdb format")
    print("====================================================================")
