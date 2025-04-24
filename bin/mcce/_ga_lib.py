"""
General Genetic Algorithm Library
"""
import os
import random
import logging
from collections import deque
from .pdbio import Residue, Conformer, Protein, Atom
from .constants import *
from ._make_connect import get_atom_by_name


class _Atom:
    """
    Internal Atom class for GA
    """
    def __init__(self):
        self.atomname = None
        self.xyz = None
        self.r_vdw = None
        self.e_vdw = None
        self.r_boundary = None
        self.charge = None
        self.connect12 = []
        self.connect13 = []

    def clone(self):
        """
        Clone this atom
        """
        new_atom = _Atom()
        new_atom.atomname = self.atomname  # atom name is immutable
        new_atom.xyz = self.xyz
        new_atom.r_vdw = self.r_vdw
        new_atom.e_vdw = self.e_vdw
        new_atom.r_boundary = self.r_boundary
        new_atom.charge = self.charge
        return new_atom

    def inherit(self, atom):
        """
        Inherit properties from an original Atom object
        """
        self.atomname = atom.atomname
        self.xyz = atom.xyz
        self.r_vdw = atom.r_vdw
        self.e_vdw = atom.e_vdw
        self.r_boundary = atom.r_boundary
        self.charge = atom.charge
        # connect12 and connect13 will be set later when creating the conformer
        # Note: connect12 and connect13 are not inherited, they will be set based on the new conformer context

    def to_original(self):
        """
        Convert this atom to the original Atom object
        """
        new_atom = Atom()
        new_atom.atomname = self.atomname
        new_atom.xyz = self.xyz
        new_atom.r_vdw = self.r_vdw
        new_atom.e_vdw = self.e_vdw
        new_atom.r_boundary = self.r_boundary
        new_atom.charge = self.charge
        new_atom.parent_conf = None
        return new_atom


class _Conformer:
    """
    Internal Conformer class for GA
    """
    def __init__(self):
        self.conftype = None
        self.atoms = []
        self.rotatables = {}
        self.parent_residue = None

    def clone(self):
        """
        Clone this conformer
        """
        new_conformer = _Conformer()
        new_conformer.conftype = self.conftype
        new_conformer.atoms = [atom.clone() for atom in self.atoms]
        new_conformer.parent_residue = self.parent_residue  # by default, a cloned conformer shares the same parent residue, update it later
        # Create a dictionary to map old atoms to new atoms
        atom_mapping = {}
        for old_atom, new_atom in zip(self.atoms, new_conformer.atoms):
            atom_mapping[old_atom] = new_atom
        # Update connect12 and connect13 for the new atoms, using the mapping for in-conformer atoms and old atoms for out-conformer atoms (backbone is shared)
        for old_atom, new_atom in atom_mapping.items():
            new_atom.connect12 = [
                atom_mapping[atom] if atom in atom_mapping else atom 
                for atom in old_atom.connect12
            ]
            new_atom.connect13 = [
                atom_mapping[atom] if atom in atom_mapping else atom 
                for atom in old_atom.connect13
            ]
        # update rotatables
        for key, value in self.rotatables.items():
            if key[0] in atom_mapping:
                atom1 = atom_mapping[key[0]]
            else:
                atom1 = key[0]
            if key[1] in atom_mapping:
                atom2 = atom_mapping[key[1]]
            else:
                atom2 = key[1]
            new_key = (atom1, atom2)
            new_conformer.rotatables[new_key] = [atom_mapping[atom] for atom in value]  # affected atoms are always in the same side chain conformer

        return new_conformer

    def inherit(self, conformer):
        """
        Inherit properties from an original Conformer object
        """
        self.conftype = conformer.conftype
        for atom in conformer.atoms:
            new_atom = _Atom()
            new_atom.inherit(atom)
            self.atoms.append(new_atom)
        self.parent_residue = conformer.parent_residue  # by default, a new conformer shares the same parent residue, in GA we will update it later
        # create a mapping old atoms to new atoms
        atom_mapping = {}
        for old_atom, new_atom in zip(conformer.atoms, self.atoms):
            atom_mapping[old_atom] = new_atom
        # update connect12 and connect13 for the new atoms
        for old_atom, new_atom in atom_mapping.items():
            new_atom.connect12 = [
                atom_mapping[atom] if atom in atom_mapping else atom 
                for atom in old_atom.connect12
            ]
            new_atom.connect13 = [
                atom_mapping[atom] if atom in atom_mapping else atom 
                for atom in old_atom.connect13
            ]
        # update rotatables
        for key, value in conformer.rotatables.items():
            if key[0] in atom_mapping:
                atom1 = atom_mapping[key[0]]
            else:
                atom1 = key[0]
            if key[1] in atom_mapping:
                atom2 = atom_mapping[key[1]]
            else:
                atom2 = key[1]
            new_key = (atom1, atom2)
            self.rotatables[new_key] = [atom_mapping[atom] for atom in value]  # affected atoms are always in the same side chain conformer

    def to_original(self):
        """
        Convert this conformer to the original Conformer objects
        """
        new_conformer = Conformer()
        new_conformer.conftype = self.conftype
        new_conformer.parent_residue = None  # set the parent residue to None for the original conformer
        for atom in self.atoms:
            new_atom = atom.to_original()
            new_atom.parent_conf = new_conformer  # set the parent conformer to the original conformer
            new_conformer.atoms.append(new_atom)            
        return new_conformer

class _Residue:
    """
    Internal Residue class for GA
    """
    def __init__(self):
        self.conformers = []
        self._flag = "GA" # flag for this residue, GA means this is a GA residue

    def clone(self):
        """
        Clone this residue
        """
        new_residue = _Residue()
        new_residue.conformers.append(self.conformers[0])  # backbone is shared
        new_residue.conformers.append(self.conformers[1].clone())  # clone the side chain conformer
        new_residue.conformers[1].parent_residue = new_residue  # set the parent residue to this residue for the new conformer
        return new_residue


class Individual:
    """
    Class for Individual in Genetic Algorithm, equivalent to a chromosome in GA
    """
    def __init__(self, pool=None):
        """
        Keep the data structure as simple as possible for performance
        """
        self.parent_pool = pool
        self.chromosome = []  # a list of flippable residues
        self.fitness = 0.0  # fitness score for this individual
        self.rank = 0  # rank for this individual (0 based)



    def create(self):
        """
        Create an individual with side chain conformers having _Atom objects instead of the original Atom objects.
        """
        for i in self.parent_pool.index_flipper:
            residue = self.parent_pool.mcce.protein.residues[i]
            selected_conformer = _Conformer()
            selected_conformer.inherit(random.choice(residue.conformers[1:]))
            selected_conformer.parent_residue = new_residue = _Residue()

            for atom in selected_conformer.atoms:
                connected_atomnames = self.parent_pool.mcce.tpl.get(
                    ("CONNECT", atom.atomname, selected_conformer.conftype), []
                ).connected
                atom.connect12 = [
                    get_atom_by_name(selected_conformer, name) or 
                    get_atom_by_name(residue.conformers[0], name)
                    for name in connected_atomnames if name
                ]
                atom.connect12 = list(filter(None, atom.connect12))
                atom.connect13 = [
                    atom2 for atom1 in atom.connect12
                    for atom2 in atom1.connect12
                    if atom2 not in atom.connect12 and atom2 != atom
                ]

            new_residue.conformers = [residue.conformers[0], selected_conformer]
            self.chromosome.append(new_residue)

    def clone(self):
        """
        Clone this individual
        """
        new_individual = Individual(pool=self.parent_pool)
        new_individual.chromosome = [residue.clone() for residue in self.chromosome]
        return new_individual
    

    def print_connect12(self):
        """
        print connect12 for this individual, debug only
        """
        for res in self.chromosome:
            for atom in res.conformers[1].atoms:
                print(atom.atomname, [a.atomname for a in atom.connect12])


    def to_mccepdb(self, fname):
        """
        Convert this individual to a mccepdb object
        """
        individual_protein = Protein()
        individual_protein.residues = [None] * len(self.parent_pool.mcce.protein.residues)

        # Restore fixed residues
        for i in self.parent_pool.index_fixed:
            individual_protein.residues[i] = self.parent_pool.mcce.protein.residues[i]

        # Restore flipper residues
        for i, residue in zip(self.parent_pool.index_flipper, self.chromosome):
            original_residue = self.parent_pool.mcce.protein.residues[i]
            individual_protein.residues[i] = residue
            residue.resname = original_residue.resname
            residue.chain = original_residue.chain
            residue.resid = original_residue.resid
            residue.sequence = original_residue.sequence
            residue.insertion = original_residue.insertion
            residue.conformers[1] = residue.conformers[1].to_original()
            residue.conformers[1].history = f"{residue.conformers[1].conftype[-2:]}G{'_'*7}"
            residue.conformers[1].parent_residue = original_residue

        individual_protein.prepend_lines = [f"# Fitness = {self.fitness:.2f}; Rank = {self.rank:d}\n"]
        individual_protein.dump(fname)
               
        
        individual_protein.prepend_lines = [f"# Fitness = {self.fitness:.2f}; Rank = {self.rank:d}\n"]
        individual_protein.dump(fname)

    def test_lineage(self):  # for debugging
        """
        Test the lineage of this , conf[0] should have parent residue as the original residue while conf[1] should have parent residue as the new residue
        """
        for res in self.chromosome:
            print(f"{res._flag} (=GA)-> {res.conformers[0].parent_residue._flag} (=Empty), {res.conformers[1].parent_residue._flag} (=GA)")

    from ._ga_forcefield import atom_embedding_depth

    def get_fitness(self):
        """
        Calculate the fitness of this individual with fast force field
        The fitness score is energy with these components:
        - vdW energy
        - elec energy
        - solvation energy, pH and Eh are NOT part of the fitness, as we are only exploring in the rotational space
        """
        # initialize the sas for atoms and side chanin conformers. *SAS relies on r_vdw to determine the SAS radius*
        # background atoms from the fixed residues
        background_atoms = []
        for i in self.parent_pool.index_fixed:
            for conf in self.parent_pool.mcce.protein.residues[i].conformers:
                background_atoms += conf.atoms
        # background atoms from the conf[0] of flipper residues
        for i in self.parent_pool.index_flipper:  # remember, conf[0] keeps lineages to the original residue
            background_atoms += self.parent_pool.mcce.protein.residues[i].conformers[0].atoms

        # get embedding depth for all atoms
        for res in self.chromosome:
            for atom in res.conformers[1].atoms:
                other_atoms = list(set(res.conformers[1].atoms) - {atom})


        # get atom list

        # loop over all atoms against all atoms to calculate the energy

            # calculate the distance

            # calculate vdw

            # calculate elec

            # compose fitness score


class Pool:
    """
    Class for Pool in Genetic Algorithm
    A pool has:
    - size: number of individuals in the pool, also equals len(population)
    - index_fixed: indices of fixed residues
    - index_flipper: indices of flipper residues
    - population: list of individuals in the pool
    """
    def __init__(self, mcce, size):
        self.mcce = mcce
        self.size = size
        self.index_fixed, self.index_flipper = self.divide_fixed_flipper()
        self.population = []
        self.pfa_queue = deque(maxlen=GA_PFA_queue)  # keep the last pfa values
        self.mcce.reset_connect()  # reset all connect12 and connect13
        self.mcce.assign_qr()  # assign charges and radii
        for i in range(size):
            individual = Individual(pool=self)
            individual.create()
            self.population.append(individual)
    
    def divide_fixed_flipper(self):
        index_fixed = []
        index_flipper = []
        for i, residue in enumerate(self.mcce.protein.residues):
            if len(residue.conformers) == 1:  # only backbone
                index_fixed.append(i)
            elif len(residue.conformers) == 2:  # one side chain
                if residue.conformers[1].rotatables:
                    index_flipper.append(i)
                else:  # one side chain but no rotatable
                    index_fixed.append(i)
            else:  # more than one side chain, flipper includes conftype
                index_flipper.append(i)
        return index_fixed, index_flipper

    def writepdb(self):
        """
        Write individuals in the pool to pdb files
        """
        # create GA_OUTPUT_FOLDER if not exists
        if not os.path.exists(GA_OUTPUT_FOLDER):
            os.makedirs(GA_OUTPUT_FOLDER)
        else:
            # remove all files in GA_OUTPUT_FOLDER
            for file in os.listdir(GA_OUTPUT_FOLDER):
                file_path = os.path.join(GA_OUTPUT_FOLDER, file)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}. Reason: {e}")
        # write individuals to pdb files
        for i, individual in enumerate(self.population[:self.size//2]):
            fname = os.path.join(GA_OUTPUT_FOLDER, f"state_{i:04d}.pdb")
            individual.to_mccepdb(fname)

    def evolve(self):
        """
        Evolve the pool
        """
        # initial fitness calculation
        for individual in self.population:
            individual.get_fitness()
        self.population.sort(key=lambda x: x.fitness)
        self.pfa = sum([individual.fitness for individual in self.population[:self.size//2]]) / (self.size//2)




def test_clone(pool):
    for individual in pool.population:
        new_individual = individual.clone()
        # Test the clone function, debug only
        # 1. check if the cloned individual has the same atoms in the residues of its chromosome
        # 2. check if the cloned atoms have correct connect12 and connect13
        # 3. check if the residues in the cloned individual have the rotatables
        if len(new_individual.chromosome) != len(individual.chromosome):
            logging.error("         The length of the chromosome is not the same after cloning.")
        for i in range(len(individual.chromosome)):
            res = individual.chromosome[i]
            new_res = new_individual.chromosome[i]
            for atom1, atom2 in zip(res.conformers[0].atoms, new_res.conformers[0].atoms):
                if atom1 != atom2:   # backbone atoms point to the same atom
                    logging.error("         The backbone atoms in the original and cloned residues are not the same.")
            for atom1, atom2 in zip(res.conformers[1].atoms, new_res.conformers[1].atoms):
                if atom1 == atom2:
                    logging.error("         The side chain atoms in the original and cloned residues point to the same atom.")
                else:
                    if atom1.atomname != atom2.atomname:
                        logging.error("         The side chain atoms in the original and cloned residues have different names.")
                    for a, b in zip(atom1.connect12, atom2.connect12):
                        if a.atomname != b.atomname:
                            logging.error("         The connect12 atoms in the original and cloned residues have different names.")
                    for a, b in zip(atom1.connect13, atom2.connect13):
                        if a.atomname != b.atomname:
                            logging.error("         The connect13 atoms in the original and cloned residues have different names.")
            zipped_keys = zip(res.conformers[1].rotatables.keys(), new_res.conformers[1].rotatables.keys())
            for key1, key2 in zipped_keys:
                if key1[0] is not None and (key1[0].atomname != key2[0].atomname):
                    logging.error("         The rotatable bond atom 1 in the original and cloned residues have different names.")
                if key1[1].atomname != key2[1].atomname:
                    logging.error(f"         The rotatable bond atom 2 in the original {key1[1].atomname} and cloned residues {key2[1].atomname} have different names.")
                affected_atoms1_atomnames = [a.atomname for a in res.conformers[1].rotatables[key1]]
                affected_atoms2_atomnames = [a.atomname for a in new_res.conformers[1].rotatables[key2]]
                
                if affected_atoms1_atomnames != affected_atoms2_atomnames:
                    logging.error("         The affected atoms in the original and cloned residues have different names.")
                    logging.error(f"         Original: {affected_atoms1_atomnames}")
                    logging.error(f"         Cloned: {affected_atoms2_atomnames}")
