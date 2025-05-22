#!/usr/bin/env python
"""
MCCE4 Tool: Compile Electrostatic Energy
Description: Compile electrostatic energy for a given protein structure (microstate)
This script does the Step 3 of the ele fitting

Step 0. Calculate the embedding score

Step 1. Prepare step2_out.pdb by ele_setup.py
delphi_ele.py will prepare a file from step2_out.pdb in which:
1. the atom radius will be set to the values used by embedding depth score calculation.
2. the atom charge will be set to +1 per conformer on one atom, so that the opp reports atom to atom ele

Step 2. Run delphi
After step2_out.pdb is ready, we will switch to MCCE4 to use step3.py run delphi
```
step3.py -s delphi --debug
```

Step 3. Compile ele from delphi and embedding score to a csv file
ele_compile.py will grab information from /tmp/pbe files and embedding score to make a csv file

Columns:
- distance between atom 1 and atom 2
- embedding score atom 1
- embedding score atom 2
- internal epsilon
- external epsilon 
- Columbs potential
- electrostatic energy calculated by delphi
"""

import logging
import argparse
import os
from mcce.geom import *

class AtomProperties:
    def __init__(self):
        self.xyz = Vector()
        self.radius = 0.0
        self.charge = 0.0
        self.embedding = 0.0

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.radius:8.3f}{self.embedding:8.3f}"

if __name__ == "__main__":
    # Set up logging    
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'

    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)    

    # Parse command line arguments
    helpmsg = """Compile embedding score and electrostatic energy calculated delphi. 
    This script needs to be called in a working directory that contains:
    1. embedding score file: statename.embedding
    2. charge file: step2_out.pdb (this file contains the charge information for each atom)
    3. electrostati energy: PBE folders under /tmp which were created by step3.py with "--debug" option

The output CSV contains columns such as distances, embedding scores, internal/external epsilon, Coulomb potential, and electrostatic energy."""
    
    parser = argparse.ArgumentParser(description=helpmsg, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("statename", nargs="?", default="", help="state name from step 2 states without file extension")
    args = parser.parse_args()

    if args.statename == "":
        logging.error("Please provide a state name from step 2 states without file extension")
        exit(1)
    
    # Load the embedding score file
    embedding_file = f"{args.statename}.embedding"
    if not os.path.isfile(embedding_file):
        logging.error(f"Embedding score file {embedding_file} not found")
        exit(1)
    with open(embedding_file, "r") as f:
        embedding_lines = f.readlines()
    embedding_lines = [line.strip() for line in embedding_lines if line.strip()]  # Remove empty lines    
    # we will store the embedding score in atom properties dictionary, indexed by atom id
    atom_properties = {}
    for line in embedding_lines:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            atomname = line[12:16]
            resname = line[17:20]
            chainid = line[21]
            resseq = line[22:26]
            atom_id = f"{atomname}_{resname}_{chainid}_{resseq}"
            atom = AtomProperties()
            atom.xyz = Vector((float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())))
            atom.radius = float(line[54:62].strip())
            atom.embedding = float(line[62:].strip())
            # print(f"{atom_id}: {atom.xyz}, {atom.radius:6.3f}, {atom.charge:6.3f}, {atom.embedding:6.3f}")
            atom_properties[atom_id] = atom    
    
    # update the atom properties with charge information from step2_out.pdb
    pdb_file = "step2_out.pdb"
    if not os.path.isfile(pdb_file):
        logging.error(f"PDB file {pdb_file} not found")
        exit(1)
    with open(pdb_file, "r") as f:
        pdb_lines = f.readlines()
    pdb_lines = [line.strip() for line in pdb_lines if line.strip()]  # Remove empty lines
    for line in pdb_lines:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            atomname = line[12:16]
            resname = line[17:20]
            chainid = line[21]
            resseq = line[22:26]
            atom_id = f"{atomname}_{resname}_{chainid}_{resseq}"
            if atom_id in atom_properties:
                atom_properties[atom_id].charge = float(line[54:62].strip())
                # print(f"{atom_id}: {atom_properties[atom_id].xyz}, {atom_properties[atom_id].radius:6.3f}, {atom_properties[atom_id].charge:6.3f}, {atom_properties[atom_id].embedding:6.3f}")