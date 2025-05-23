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

# Constants
k_coulomb = 332.06371  # Coulomb's constant in kcal/(mol*Å*e^2)


class AtomProperties:
    def __init__(self):
        self.confid = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.charge = 0.0
        self.embedding = 0.0

    def __repr__(self):
        return f"{self.line[:30]}{self.xyz.x:8.3f}{self.xyz.y:8.3f}{self.xyz.z:8.3f}{self.radius:8.3f}{self.embedding:8.3f}"

def get_embedding_score(fname):
    """Create atom_properties dictionary and initialize it with atom embedding score."""
    atom_properties = {}

    # Load the embedding score file
    if not os.path.isfile(fname):
        logging.error(f"Embedding score file {fname} not found")
        exit(1)
    embedding_lines = open(fname).readlines()
    embedding_lines = [line.strip() for line in embedding_lines if line.strip()]  # Remove empty lines    
    # we will store the embedding score in atom properties dictionary, indexed by atom id
    for line in embedding_lines:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            atomname = line[12:16]
            resname = line[17:20]
            chainid = line[21]
            resseq = line[22:26]
            atom_id = (atomname, resname, chainid, resseq)
            atom = AtomProperties()
            atom.xyz = Vector((float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())))
            atom.radius = float(line[54:62].strip())
            atom.embedding = float(line[62:].strip())
            # print(f"{atom_id}: {atom.xyz}, {atom.radius:6.3f}, {atom.charge:6.3f}, {atom.embedding:6.3f}")
            atom_properties[atom_id] = atom    
    
    return atom_properties

def update_charge(atom_properties):
    """Update the atom properties with charge information from step2_out.pdb."""
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
            atom_id = (atomname, resname, chainid, resseq)
            confid = resname + line[80:82] + line[21:30]
            if atom_id in atom_properties:
                atom_properties[atom_id].charge = float(line[62:74].strip())
                atom_properties[atom_id].confid = confid


def get_electrostatic_energy(atom_properties):
    """
    Get electrostatic energy from PBE files.
    Electrostatic energy is calculated by delphi and stored in /tmp/pbe files.
    Coulomb potential is calculated by delphi and stored in /tmp/pbe files.
    PBE files will be loaded once and stored in a dictionary.
    What goes into the csv file is the matrix of atoms from atom_proterties if the atom charge is not 0.  
    """
    # Get the PBE file name
    cwd = os.getcwd().replace("/", ".").strip(".")
    prefix = f"/tmp/pbe_{cwd}."

    for atom_id, atom in atom_properties.items():
        if atom.charge == 0:
            continue
        pbe_file = prefix+atom.confid
        if not os.path.isdir(pbe_file):
            logging.error(f"PBE folder {pbe_file} not found")
            exit(1)




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
    
    logging.info(f"Loading embedding score from {args.statename}.embedding ...")
    atom_properties = get_embedding_score(f"{args.statename}.embedding")
    logging.info("Updating charge and conformer ID for atoms using step2_out.pdb ...")
    update_charge(atom_properties)
    logging.info("Obtaining Calculating electrostatic energy from pbe folders under /tmp ...")
    electrostatic_energy = get_electrostatic_energy(atom_properties)
