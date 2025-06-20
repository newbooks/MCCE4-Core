#!/usr/bin/env python
"""
MCCE4 Tool: Compile Electrostatic Energy
Description: Compile electrostatic energy for a given protein structure (microstate)
This script does the Step 3 of the ele fitting

Required files:
- step2_out.pdb: PDB file with atom properties (charge, radius, etc.)
- statename.density: Local density file for the state name
- energies/step3_out.opp: Opp file with electrostatic energy calculated by delphi

Output file:
Columns:
- Conf1: Conformer ID for atom 1
- Conf2: Conformer ID for atom 2
- Distance: Distance between two atoms in Angstroms
- Radius1: radius of atom 1
- Radius2: radius of atom 2
- Density1_Near: Near density score for atom 1 
- Density2_Near: Near density score for atom 2
- Density1_Mid: Far density score for atom 1
- Density2_Mid: Far density score for atom 2
- Density1_Far: Far density score for atom 1
- Density2_Far: Far density score for atom 2
- CoulombPotential: Coulomb potential in kcal/mol
- PBPotential: electrostatic energy calculated by delphi in kcal/mol
"""

import logging
import argparse
import os
from mcce.geom import *


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

    logging.info(f"Loading local density from {args.statename}.density ...")
    update_density_score(atoms, f"{args.statename}.density")

    # print the atom properties
    # for atom_id, atom in atoms.items():
    #     logging.info(f"Atom ID: {atom_id}, Confid: {atom.confid}, XYZ: {atom.xyz}, Radius: {atom.radius}, Charge: {atom.charge}, Embedding: {atom.embedding}")


    logging.info("Obtaining electrostatic energy from opp ...")
    pairwise_ele = get_electrostatic_energy(atoms)

    logging.info("Compiling results into CSV file ...")
    output_file = f"{args.statename}_compiled.csv"
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