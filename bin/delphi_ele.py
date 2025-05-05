#!/usr/bin/env python
"""
Program Name: delphi_ele
Description: Set up a single charge on each conformer and calculate the electrostatic energy
Output:
    energies/*.opp  for electrostatic energy
"""

import logging
import argparse
from mcce.geom import *

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

# Make sure this is consistent with the one in embedding_score.py
ATOM_RADII = {
    ' H': 1.2,
    ' C': 1.7,
    ' N': 1.55,
    ' O': 1.52,
    ' S': 1.8,
    ' P': 1.8,
    # Add more atom types and their radii as needed
}
ATOM_RADIUS_UNKNOWN = 1.5  # Default radius for unknown atom types


def parse_arguments():
    helpmsg = "Set up step2_out.pdb for calculating ele by delphi"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", nargs='?', default="", help="PDB file to be processed")
    
    return parser.parse_args()


class Atom:
    def __init__(self):
        self.line = ""  # Original line from PDB file
        self.element = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.charge = 0.0


def read_pdb(pdb_file):
    """
    Read a PDB file and extract atom information.
    """
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                atom = Atom()
                atom.line = line.strip()
                atomname = line[12:16]
                if len(atomname.strip()) == 4 and atomname[0] == "H":
                    atom.element = " H"
                else:
                    atom.element = atomname[:2]
                atom.xyz = Vector((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                atom.radius = ATOM_RADII.get(atom.element, ATOM_RADIUS_UNKNOWN)

                atom.xyz = Vector(float(line[30:38]), float(line[38:46]), float(line[46:54]))
                atom.radius = ATOM_RADII.get(atom.element, ATOM_RADIUS_UNKNOWN)
                atoms.append(atom)
    return atoms


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)
    args = parse_arguments()
    if args.pdb_file:
        logging.info(f"Processing PDB file: {args.pdb_file}")
    else:
        logging.error("No PDB file provided.")
        exit(1)

    atoms = read_pdb(args.pdb_file)
