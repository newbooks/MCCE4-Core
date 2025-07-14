#!/usr/bin/env python3
"""
All in one training script for fast force field by Machine Learning.
"""

import logging
import argparse
import os
import shutil
import subprocess
import tempfile
import random
import time
from mcce.geometry import Vector


# Atom radii are for step2_out.pdb
ATOM_RADII = {
    " H": 1.20,  # Hydrogen
    " C": 1.70,  # Carbon
    " N": 1.55,  # Nitrogen
    " O": 1.52,  # Oxygen
    " S": 1.80,  # Sulfur
    " P": 1.80,  # Phosphorus
}
ATOM_RADIUS_UNKNOWN = 1.80  # Default radius for unknown elements

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
    Assign charges to atoms of the protein. +1 charge is assigned to a random atom of each side chain conformer as long as the atom.
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

class AtomProperties:
    __slots__ = ("confid", "xyz", "radius", "charge", "density_near", "density_mid", "density_far", "d2surface")

    def __init__(self):
        self.confid = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.charge = 0.0
        self.density_near = 0
        self.density_mid = 0
        self.density_far = 0
        self.d2surface = 0.0 

    def __repr__(self):
        return (f"{self.confid} {self.xyz.x:8.3f} {self.xyz.y:8.3f} {self.xyz.z:8.3f} "
                f"{self.radius:8.3f} {self.charge:8.3f} {self.density_near:6d} {self.density_mid:6d} "
                f"{self.density_far:6d} {self.d2surface:6.3f}")

def load_atoms():
    """
    Load atom properties from step2_out.pdb.
    Only charged atoms are considered.
    """
    pdb_file = "step2_out.pdb"
    if not os.path.isfile(pdb_file):
        logging.error(f"PDB file {pdb_file} not found")
        exit(1)
    atoms = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atom = AtomProperties()
                atomname = line[12:16]
                resname = line[17:20]
                chainid = line[21]
                resseq = line[22:26]
                atom_id = (atomname, resname, chainid, resseq)
                confid = resname + line[80:82] + line[21:30]
                atom.xyz = Vector((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                atom.radius = float(line[54:62])
                atom.charge = float(line[62:74])
                if abs(atom.charge) > 0.001:  # Only consider charged atoms
                    atom.confid = confid
                    atoms[atom_id] = atom
    return atoms



if __name__ == "__main__":
    # Set up command line arguments
    helpmsg = "Train a fast force field model from a pdb folder. This script attempts to reproduce PB solver delphi potential using machine learning techniques."
    parser = argparse.ArgumentParser(description=helpmsg)
    parser.add_argument('pdb_folder', type=str, help='Path to the pdb folder. The folder should contain proteins at different size and shape.')
    args = parser.parse_args()

    # Set up logging
    logging_format = "%(asctime)s %(levelname)s: %(message)s"
    logging_datefmt='%Y-%m-%d %H:%M:%S'
    logging.basicConfig(format=logging_format, datefmt=logging_datefmt, level=logging.INFO)

    # Check if the provided path is a directory and contains PDB files
    if not os.path.isdir(args.pdb_folder):
        logging.error(f"The provided path '{args.pdb_folder}' is not a valid directory.")
        exit(1)
    if not any(file.endswith('.pdb') for file in os.listdir(args.pdb_folder)):
        logging.error(f"The directory '{args.pdb_folder}' does not contain any PDB files.")
        exit(1)
    
    # get the absolute path to pdb files in the folder
    args.pdb_folder = os.path.abspath(args.pdb_folder)
    pdbs = [os.path.join(args.pdb_folder, file) for file in os.listdir(args.pdb_folder) if file.endswith('.pdb')]

    # Work under a temporary directory
    logging.info(f"Using temporary directory for training: {args.pdb_folder}")
    with tempfile.TemporaryDirectory(prefix="aiofff_training_") as temp_dir:
        logging.info(f"Created temporary directory at {temp_dir}")
        # go to the temporary directory
        os.chdir(temp_dir)
        for pdb in pdbs[:1]: # limit to the first pdb file for testing
            base_name = os.path.basename(pdb)
            logging.info(f"Processing PDB file: {base_name}")
            
            logging.info(f"Running step1.py")
            shutil.copy(pdb, temp_dir)
            result = subprocess.run([f"step1.py", f"{base_name}"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"step1.py failed on {base_name} with exit code {result.returncode}")
                exit(result.returncode)

            logging.info(f"Running step2.py")
            result = subprocess.run([f"step2.py", "--writepdb"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"step2.py failed on {base_name} with exit code {result.returncode}")
                exit(result.returncode)
            shutil.copy("ga_output/state_0001.pdb", "microstate.pdb")
            
            logging.info(f"Running local_density.py")
            result = subprocess.run([f"local_density.py", "microstate.pdb"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"local_density.py failed on microstate.pdb with exit code {result.returncode}")
                exit(result.returncode)

            logging.info(f"Setting charge and radius and write step2_out.pdb")
            atoms = read_pdb("microstate.pdb")
            assign_charges(atoms)
            write_pdb(atoms, "step2_out.pdb")

            logging.info(f"Running step3.py")
            result = subprocess.run([f"step3.py", "-s", "delphi", "-p", "3"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                logging.error(f"step3.py failed with exit code {result.returncode}")
                exit(result.returncode)

            logging.info("Compile density and electrostatics energy to a csv file.")
            atoms = load_atoms()