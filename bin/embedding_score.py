#!/usr/bin/env python

"""
MCCE4 Tool: Standalone Program to Calculate Atom Embedding Score
"""

import logging
import argparse
from mcce.geom import *

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

# Constants
PROBE_RAD = 1.4  # Probe radius in Angstrom
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
    helpmsg = "Calculate atom embedding score for a given protein structure"
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument("pdb_file", nargs='+', default=[], help="PDB file(s) to be processed")
    
    return parser.parse_args()


class Atom:
    def __init__(self):
        self.line = ""  # Original line from PDB file
        self.element = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.embedding = 0.0

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.radius:8.3f}{self.embedding:8.3f}"

def load_atoms_from_pdb(pdb_file):
    lines = open(pdb_file).readlines()
    atoms = []
    
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom = Atom()
            atom.line = line.strip()
            atomname = line[12:16]
            if len(atomname.strip()) == 4 and atomname[0] == "H":
                atom.element = " H"
            else:
                atom.element = atomname[:2]
            atom.xyz = Vector((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            atom.radius = ATOM_RADII.get(atom.element, ATOM_RADIUS_UNKNOWN)
            atoms.append(atom)
   
    return atoms



if __name__ == "__main__":
    args = parse_arguments()
    for pdb in args.pdb_file:
        print(f"Processing {pdb}...")            
        atoms = load_atoms_from_pdb(pdb)
    for atom in atoms:
        print(atom)