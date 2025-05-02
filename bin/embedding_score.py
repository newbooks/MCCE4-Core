#!/usr/bin/env python

"""
MCCE4 Tool: Standalone Program to Calculate Atom Embedding Score
"""

import logging
import argparse
from .mcce.geom import *

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

    def __repr__(self):
        return f"Atom(type={self.atom_type}, coords={self.coordinates}, radius={self.radius})"

def load_atoms_from_pdb(pdb_file):
    pass


if __name__ == "__main__":
    args = parse_arguments()
    for pdb in args.pdb_file:
        print(f"Processing {pdb}...")            
        atoms = load_atoms_from_pdb(pdb)