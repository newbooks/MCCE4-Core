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
import random
import time
# Seed the random number generator with the current time
random.seed(time.time())

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
                atom.radius = ATOM_RADII.get(atom.element, ATOM_RADIUS_UNKNOWN)

                atoms.append(atom)
    return atoms

def assign_charges(atoms):
    """
    Assign charges to atoms of the protein. +1 charge is assigned to a random atom of each side chain conformer.
    The reason a random atom is chosen is that more atoms with different embedding depth will be sampled.
    """

    sidechain_atoms = {}
    for atom in atoms:
        # colelct all atoms of the same side chain
        sidechain_id = atom.line[17:30]
        if sidechain_id[-3:] != "000":
            if sidechain_id not in sidechain_atoms:
                sidechain_atoms[sidechain_id] = []
            sidechain_atoms[sidechain_id].append(atom)
    for sidechain_id, sidechain_atoms in sidechain_atoms.items():
        random_atom = random.choice(sidechain_atoms)
        random_atom.charge = 1.0
        

def write_pdb(atoms, output_file):
    """
    Write the modified atoms to a new PDB file.
    """
    with open(output_file, 'w') as f:
        for atom in atoms:
            line = atom.line[:54] + f"{atom.radius:8.3f}" + f"{atom.charge:12.3f}" + atom.line[74:]
            f.write(line + '\n')

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)
    args = parse_arguments()
    if args.pdb_file:
        logging.info(f"Processing PDB file: {args.pdb_file}")
    else:
        logging.error("No PDB file provided.")
        exit(1)

    atoms = read_pdb(args.pdb_file)
    assign_charges(atoms)
    write_pdb(atoms, "step2_out.pdb")
    logging.info(f"MCCE PDB file written to step2_out.pdb for step3.py")