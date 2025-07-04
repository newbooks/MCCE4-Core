#!/usr/bin/env python

"""
MCCE4 Step 1: Read PDB file and generate a protein object
"""

import logging
import argparse
from mcce.pdbio import *
from mcce.main import MCCE

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
    parser.add_argument("--no_center", default=False, action="store_true", help="This option prevents moving the weighted center of protein to the origin.")
    parser.add_argument("--cofactors", metavar="cofactor", nargs="+", default=[], help="Load additional cofactors that can be stripped off, in quotes if space in name")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    if args.cofactors:
        LOOSE_COFACTORS += args.cofactors

    logging.info("Step 1: Read PDB file and generate a mcce protein object")
    
    # Get prm
    prm = Runprm()                  # Create a Runprm instance and load the default runprm file
    prm.update_by_cmd(args)         # Update parameters using command line arguments
    prm.update_by_files(args.r)     # Update parameters using additional runprm files
    prm.dump(comment="Step 1 uses these runprm parameters") # Save the parameters to run.prm.record

    # Get tpl
    tpl = Tpl()                     # Create a Tpl instance
    tpl.load_ftpl_folder(prm._FTPL_FOLDER.value) # Load the ftpl folder specified in runprm
    if os.path.isdir(USER_PARAM):
        tpl.load_ftpl_folder(USER_PARAM)  # Load user defined ftpl files
    tpl.dump(comment="Step 1 uses these tpl parameters") # Save the parameters to tpl.dat.record

    # Get pdb
    # Attempt to load the pdb file to a Pdb object.
    # The Pdb object has these attributes that are relavent to the next part of step 1:
    #   - pdb.mcce_ready: a boolean value indicating whether the pdb file is ready for mcce
    #   - pdb.message: a string message indicating the why the pdb file is not mcce ready
    #   - pdb.atoms: a list of Atom objects
    pdb = Pdb(prm.INPDB.value)  # Create a Pdb instance and load the pdb file to atoms
    if not pdb.mcce_ready:
        logging.error(f"{pdb.message}")
        exit(1)

    if args.no_hoh:
        logging.info("   Removing water molecules")
        pdb.remove_hoh()  # Remove water molecules

    pdb.identify_ligands(tpl)   # Identify ligands in the pdb file
    # pdb.dump_pdb("debug.pdb")       # Save the pdb to a string
    pdb.rename(prm._RENAME_RULES.value)  # Rename atoms according to the rules in runprm


    # Convert the pdb to a Protein object
    protein = pdb.convert_to_protein(tpl)  # Convert the pdb to a Protein object
    if prm.TERMINALS.value.lower() == "t":
        logging.info("   Making terminal residues")
        protein.make_ter_residues()  # Add terminal residues if necessary
    
    # Now the protein object has the right residue names, we will 
    # 1. scan protein residues to find unknown cofactors, if found, create new_ftpl for them and ammend tpl database
    # 2. move backbone atoms to conformer[0] or create new_ftpl for unknown cofactors
    # 3. split altloc atoms to conformers
    # 4. assign conformer types to conformers
    protein.new_ftpl(tpl)  # Create a new_ftpl for unknown cofactors, if any
    protein.split_backbone(tpl)  # Split backbone atoms to conformer[0]
    protein.split_altloc()  # Split altloc atoms to conformers
    # for res in protein.residues:
    #     for conf in res.conformers:
    #         print(conf.history)
    protein.assign_conftype(tpl)  # Assign conformer types to conformers

    # By this point, we gave
    # prm - Runprm object about how to run MCCE
    # tpl - Tpl object about the force field and structure paramters for pronteins
    # protein - Protein object that represents the actual input protein
    # From here, we have all three components to make a MCCE object and run MCCE.
    # Here is a test run, create a MCCE object and calculate solvent accessible surface area of the protein
    mcce = MCCE(prm=prm, tpl=tpl, protein=protein)

    while (n_stripped := mcce.strip_cofactors()) > 0:
        logging.info(f"   {n_stripped} cofactors stripped in this round")

    # calculate SAS for all residues in the protein final structure and write to acc.res and acc.atm
    mcce.sas_protein()

    if not args.no_center:
        protein.center_to_origin()

    # create head1.lst
    mcce.make_head1(STEP1_HEAD)


    protein.dump(STEP1_OUT)  # Save the protein to a pdb file
    
    # print summary
    
    print("====================================================================")
    print("Step 1 completed. Here is run summary:")
    print("--------------------------------------------------------------------")    
    print(f"  {STEP1_OUT:<16s}: Protein object in mccepdb file, used in step 2")
    print(f"  {STEP1_HEAD:<16s}: head1.lst for optional step 2 instructions")
    print(f"  {ACC_RES:<16s}: SAS for all residues")
    print(f"  {ACC_ATOM:<16s}: SAS for all atoms")
    print("====================================================================")
    