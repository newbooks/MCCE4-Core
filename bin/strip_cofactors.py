#!/usr/bin/env python

"""
Strip exposed cofactors from a pdb file
"""

import argparse
import logging
from mcce.pdbio import *
from mcce.constants import LOOSE_COFACTORS
from mcce._sas import sas_atoms
import os
import mcce._sas

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

def strip_cofactors(pdb):
    # group atoms by residue
    residues = defaultdict(list)
    for atom in pdb.atoms:
        residues[atom.residue_id()].append(atom)

    # calculate SAS for each residue
    exclude_cofactors = set()
    all_atoms = set(pdb.atoms)
    for res, atoms in residues.items():
        if res[0] in loose_cofactors:
            background = list(all_atoms - set(atoms))
            res_sas = sas_atoms(atoms, background)
            res_sas_exposed = sas_atoms(atoms, [])
            res_sas_percentage = res_sas / res_sas_exposed
            if res_sas_percentage > float(args.s):
                exclude_cofactors.add(res)
    
    # remove the cofactors from the pdb atoms list
    new_atoms = [atom for atom in pdb.atoms if atom.residue_id() not in exclude_cofactors]
    pdb.atoms = new_atoms
    return len(exclude_cofactors)


if __name__ == "__main__":
    helpmsg = "Strip off cofactors from a pdb file if the SAS is larger than the cutoff"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", metavar="pdb", nargs=1, default=[], help="PDB file")
    parser.add_argument("-s", metavar="cutoff", type=float, default=0.05, help="SAS cutoff, default is 0.05")
    parser.add_argument("--probe", metavar="probe", type=float, default=1.4, help="PDB file")
    parser.add_argument("--cofactors", metavar="cofactor", nargs="+", default=[], help="Load additional cofactors that can be stripped off, in quotes if space in name")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    pdb_file = args.pdb_file[0]
    mcce._sas.probe_global = args.probe  # update global probe radius with the input probe radius
    loose_cofactors = LOOSE_COFACTORS + args.cofactors

    # The standalone script runs on Pdb object, so the output will be a similar PDB format file rather than MCCE PDB format
    pdb = Pdb(pdb_file)
    if not pdb.mcce_ready:
        logging.error(f"{pdb.message}")
        exit(1)

    # update global probe radius with the input probe radius
    probe_global = args.probe

    # log parameters used
    logging.info(f"Probe radius: {probe_global:.2f}")
    logging.info(f"SAS cutoff: {args.s:.2f}")
    logging.info(f"Loose cofactors: {loose_cofactors}")

    # assign radius to atoms
    for atom in pdb.atoms:
        if atom.r_vdw < 0.0001:
            atom.r_vdw = mcce._sas.radius.get(atom.element, UNASSIGEDN_RAD)
            logging.debug(f"Warning: atom \"{atom.element}\" has no VDW radius, using default {atom.r_vdw}")

    while (n_stripped := strip_cofactors(pdb)) > 0:
        logging.info(f"Removed {n_stripped} cofactors this round")

    pdb_file_base = os.path.splitext(pdb_file)[0]
    pdb.dump_pdb(f"{pdb_file_base}.stripped.pdb")
    logging.info(f"Stripped pdb file saved to {pdb_file_base}.stripped.pdb")