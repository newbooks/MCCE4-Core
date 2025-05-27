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
ele_compile.py will grab information from opp files under energies folder and embedding score to make a csv file

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

def update_embedding_score(atoms, fname):
    """
    Update embedding scores for atoms.
    """
    
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
            if atom_id in atoms:
                atoms[atom_id].embedding = float(line[62:].strip())            
    
def load_atoms():
    """
    Load atom properties from step2_out.pdb.
        - confid
        - xyz
        - radius
        - charge
    We only care about the charge atoms. That is, one atom from one side chain. This way, the geometry environment of atoms is distinct.
    """
    pdb_file = "step2_out.pdb"
    if not os.path.isfile(pdb_file):
        logging.error(f"PDB file {pdb_file} not found")
        exit(1)
    pdb_lines = open(pdb_file).readlines()
    pdb_lines = [line.strip() for line in pdb_lines if line.strip()]  # Remove empty lines
    atoms = {}
    for line in pdb_lines:
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            atom = AtomProperties()
            atomname = line[12:16]
            resname = line[17:20]
            chainid = line[21]
            resseq = line[22:26]
            atom_id = (atomname, resname, chainid, resseq)
            confid = resname + line[80:82] + line[21:30]
            atom.xyz = Vector((float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())))
            atom.radius = float(line[54:62].strip())
            atom.charge = float(line[62:74].strip())
            atom.confid = confid
            if abs(atom.charge) > 0.1:
                atoms[atom_id] = atom

    return atoms


def get_electrostatic_energy(atoms):
    """
    Get electrostatic energy from PBE files.
    Electrostatic energy is calculated by delphi and stored in /tmp/pbe files.
    Coulomb potential is calculated by delphi and stored in /tmp/pbe files.
    PBE files will be loaded once and stored in a dictionary.
    What goes into the csv file is the matrix of atoms from atom_proterties if the atom charge is not 0.  
    """

    # reverse search for the atom properties
    conformers = {}
    for atom_id, atom in atoms.items():
        if atom.confid not in conformers:
            conformers[atom.confid] = atom_id
        else:
            logging.warning(f"Confid {atom.confid} already exists, skipping {atom_id}")

    # Get the PBE file name
    confids = list(conformers.keys())
    for iconf in range(len(confids)-1):
        conf1 = confids[iconf]
        fname = "energies/"+ conf1 + ".opp"
        if not os.path.isfile(fname):
            logging.error(f"Opp file {fname} not found")
            exit(1)
        opp_lines = open(fname).readlines()



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
    
    logging.info("Loading atoms from step2_out.pdb ...")
    atoms = load_atoms()

    logging.info(f"Loading embedding score from {args.statename}.embedding ...")
    update_embedding_score(atoms, f"{args.statename}.embedding")

    # print the atom properties
    # for atom_id, atom in atoms.items():
    #     logging.info(f"Atom ID: {atom_id}, Confid: {atom.confid}, XYZ: {atom.xyz}, Radius: {atom.radius}, Charge: {atom.charge}, Embedding: {atom.embedding}")


    logging.info("Obtaining electrostatic energy from opp ...")
    electrostatic_energy = get_electrostatic_energy(atoms)
