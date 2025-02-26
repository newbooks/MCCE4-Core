#!/usr/bin/env python

"""
Split structure if multiple altLoc found on backbone atoms
"""

import os
import argparse
import logging
from mcce.pdbio import *
from mcce.constants import *

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

if __name__ == "__main__":
    helpmsg = "Split structure if multiple altLoc found on backbone atoms"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_files", metavar="pdbs", nargs='+', default=[], help="PDB files")
    parser.add_argument("--debug", default=False, action="store_true", help="Print debug information")
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logging_format_debug, datefmt=logging_datefmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)

    logging.info("Split structure if multiple altLoc found on backbone atoms")

    for pdb_file in args.pdb_files:
        if os.path.isfile(pdb_file):
            atoms = []
            pdblines = open(pdb_file, "r").readlines()
            for line in pdblines:
                if line.startswith("ATOM  ") or line.startswith("HETATM"):
                    atom = Atom()
                    atom.load_pdbline(line)
                    atoms.append(atom)
            
            altloc = set()
            for atom in atoms:
                if atom.atomname in BACKBONE_ATOMS and atom.resname in RESIDUE_NAMES and atom.altloc != " ":
                    altloc.add(atom.altloc)
            
            if len(altloc) > 1:
                logging.info(f"Splitting {pdb_file}")
                altloc = sorted(list(altloc))
                fnames = []
                pdb_file_base = os.path.splitext(pdb_file)[0]
                for a in altloc:
                    with open(f"{pdb_file_base}.altloc{a}.pdb", "w") as f:
                        fnames.append(f"{pdb_file_base}.altloc{a}.pdb")
                        for line in pdblines:
                            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                                atom = Atom()
                                atom.load_pdbline(line)
                                if atom.altloc == " " or atom.altloc == a:  # Keep the atom if altloc is blank or matches the current altloc
                                    f.write(line)
                            else:
                                f.write(line)
                logging.info(f"Generated {len(altloc)} structures: {', '.join(fnames)}")
            else:
                logging.info(f"{pdb_file} has no multiple altLoc")