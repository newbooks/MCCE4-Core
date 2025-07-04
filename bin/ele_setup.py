#!/usr/bin/env python
"""
Program Name: ele_setup.py
Description: Set up a single charge on each conformer and calculate the electrostatic energy
Output:
    CSV file with columns:
    Distance, DensityAverage_Near, DensityAverage_Mid, DensityAverage_Far, P
"""

import logging
import argparse
import sys
import os
import subprocess
from mcce.geom import *
import random
import time
# Seed the random number generator with the current time
random.seed(time.time())

logging_format = "%(asctime)s %(levelname)s: %(message)s"
logging_format_debug = "%(asctime)s %(levelname)s [%(module)s]: %(message)s"
logging_datefmt='%Y-%m-%d %H:%M:%S'

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

N_THREADS = 3  # Number of threads for step3.py, can be set to a higher value if needed


def parse_arguments():
    helpmsg = "Set up and run a microstate pdb to make a csv file for ele training and validation."
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


class AtomProperties:
    __slots__ = ("confid", "xyz", "radius", "charge", "density_near", "density_mid", "density_far")

    def __init__(self):
        self.confid = ""
        self.xyz = Vector()
        self.radius = 0.0
        self.charge = 0.0
        self.density_near = 0
        self.density_mid = 0
        self.density_far = 0

    def __repr__(self):
        return (f"{self.confid} {self.xyz.x:8.3f} {self.xyz.y:8.3f} {self.xyz.z:8.3f} "
                f"{self.radius:8.3f} {self.charge:8.3f} {self.density_mid:8d} "
                f"{self.density_far:8d}")


def update_density_score(atoms, fname):
    """
    Update local density scores for atoms.
    Local density is calculated by step2_out.pdb and stored in a file with the same name as the state name.
    The file contains lines like:
    ATOM  1  N   ALA A   1      11.104  12.345  13.456  1.00  0.00           N
    The last columns are the local density scores.
    """
    if not os.path.isfile(fname):
        logging.error(f"Local density file {fname} not found")
        exit(1)
    with open(fname) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                atom_id = (line[12:16], line[17:20], line[21], line[22:26])
                local_density = [int(x) for x in line[54:].strip().split()]
                if atom_id in atoms:
                    # Expecting three density values at the end of the line
                    atoms[atom_id].density_near = local_density[0]
                    atoms[atom_id].density_mid = local_density[1]
                    atoms[atom_id].density_far = local_density[2]


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


def get_electrostatic_energy(atoms):
    """
    Get electrostatic energy from opp files.
    Returns a dictionary mapping (atom_id1, atom_id2) to electrostatic energy.
    """
    # Map confid to atom_id for quick lookup
    conformers = {}
    for atom_id, atom in atoms.items():
        if atom.confid not in conformers:
            conformers[atom.confid] = atom_id
        else:
            logging.warning(f"Confid {atom.confid} already exists, skipping {atom_id}")

    pairwise_ele = {}
    opp_dir = "energies"
    for conf1, atom_id1 in conformers.items():
        fname = os.path.join(opp_dir, f"{conf1}.opp")
        if not os.path.isfile(fname):
            logging.error(f"Opp file {fname} not found")
            exit(1)
        with open(fname) as opp_file:
            for line in opp_file:
                fields = line.strip().split()
                if len(fields) < 6:
                    continue
                conf2 = fields[1]
                atom_id2 = conformers.get(conf2)
                if atom_id2 is not None:
                    try:
                        ele = float(fields[5])
                        reverse_key = (atom_id2, atom_id1)
                        if pairwise_ele.get(reverse_key) is None:  # the other direction is not recorded
                            pairwise_ele[(atom_id1, atom_id2)] = ele
                        else:
                            pairwise_ele[reverse_key] = (pairwise_ele[reverse_key] + ele) / 2  # Average the energy if both directions are recorded
                    except ValueError:
                        continue
    return pairwise_ele



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format=logging_format, datefmt=logging_datefmt)
    args = parse_arguments()
    if args.pdb_file:
        logging.info(f"Processing PDB file: {args.pdb_file}")
    else:
        logging.error("No PDB file provided.")
        exit(1)

    # Run local density calculation
    logging.info("Running local_density.py to calculate local density ...")
    result = subprocess.run(["local_density.py", args.pdb_file], capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"local_density.py failed with exit code {result.returncode}")
        logging.error(f"stderr: {result.stderr.strip()}")
        sys.exit(1)

    # Read the PDB file and assign charges and radii
    atoms = read_pdb(args.pdb_file)
    assign_charges(atoms)
    write_pdb(atoms, "step2_out.pdb")
    # Log the completion of the PDB file writing
    logging.info("Completed writing step2_out.pdb and ready to run step3.py.")


    # Run step3.py to calculate electrostatic energy
    result = subprocess.run(["step3.py", "-s", "delphi", "-p", f"{N_THREADS}"], capture_output=True, text=True)
    if result.returncode != 0:
        logging.error(f"step3.py failed with exit code {result.returncode}")
        logging.error(f"stderr: {result.stderr.strip()}")
        sys.exit(1)
    else:
        logging.info("Completed step3.py successfully.")
    

    # Compile the results
    logging.info("Loading atoms from step2_out.pdb ...")
    atoms = load_atoms()

    pdb_base = os.path.splitext(args.pdb_file)[0]
    logging.info(f"Loading local density from {pdb_base}.density ...")
    update_density_score(atoms, f"{pdb_base}.density")

    logging.info("Obtaining electrostatic energy from opp ...")
    pairwise_ele = get_electrostatic_energy(atoms)

    logging.info("Compiling results into CSV file ...")
    output_file = f"{pdb_base}_compiled.csv"
    with open(output_file, 'w') as f:
        # We only need to write these features to the CSV file
        # Distance, DensityAverage_Near, DensityAverage_Mid, DensityAverage_Far, PBPotential
        f.write("Distance,DensityAverage_Near,DensityAverage_Mid,DensityAverage_Far,PBPotential\n")
        for (atom_id1, atom_id2), ele in pairwise_ele.items():
            atom1, atom2 = atoms[atom_id1], atoms[atom_id2]
            distance = atom1.xyz.distance(atom2.xyz)
            densityaverage_near = (atom1.density_near + atom2.density_near) / 2
            densityaverage_mid = (atom1.density_mid + atom2.density_mid) / 2
            densityaverage_far = (atom1.density_far + atom2.density_far) / 2
            # write the row to the CSV file
            f.write(f"{distance:.3f},{densityaverage_near:.3f},{densityaverage_mid:.3f},{densityaverage_far:.3f},{ele:.3f}\n")
    logging.info(f"Results compiled into {output_file}. Energy unit is kcal/mol.")