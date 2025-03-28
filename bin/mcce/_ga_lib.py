"""
General Genetic Algorithm Library
"""
import os
import random
from collections import deque
from .pdbio import Residue, Protein, Atom
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
        self.crg = None
        self.connect12 = []
        self.connect13 = []

    def clone(self):
        """
        Clone this atom
        """
        new_atom = _Atom()
        new_atom.xyz = self.xyz
        new_atom.r_vdw = self.r_vdw
        new_atom.e_vdw = self.e_vdw
        new_atom.r_boundary = self.r_boundary
        new_atom.crg = self.crg
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
        new_residue.conformers[0] = self.conformers[0]  # backbone is shared
        new_residue.conformers[1] = self.conformers[1].clone()  # clone the side chain conformer
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
        Create an individual
        """
        for i in self.parent_pool.index_flipper:
            residue = self.parent_pool.mcce.protein.residues[i]
            selected_conformer = random.choice(residue.conformers[1:]).clone()
            
            new_residue = _Residue()
            selected_conformer.parent_residue = new_residue

            for atom in selected_conformer.atoms:
                key = ("CONNECT", atom.atomname, selected_conformer.conftype)
                if key in self.parent_pool.mcce.tpl:
                    connected_atomnames = self.parent_pool.mcce.tpl[key].connected
                else:
                    connected_atomnames = []
                atom.connect12 = [
                    get_atom_by_name(selected_conformer, name) or get_atom_by_name(residue.conformers[0], name)
                    for name in connected_atomnames if name
                ]
                atom.connect12 = [a for a in atom.connect12 if a]  # filter out None

                atom.connect13 = [
                    connected_atom2
                    for connected_atom in atom.connect12
                    for connected_atom2 in connected_atom.connect12
                    if connected_atom2 not in atom.connect12 and connected_atom2 != atom
                ]

            new_residue.conformers = [residue.conformers[0], selected_conformer]
            self.chromosome.append(new_residue)

    def clone(self):
        """
        Clone this individual
        """
        new_individual = Individual(pool=self.parent_pool)
        #new_individual.chromosome = [residue.clone() for residue in self.chromosome]
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
        individual_protein.residues = len(self.parent_pool.mcce.protein.residues) * [None]
        # put back fixed residues
        for ires, residue in zip(self.parent_pool.index_fixed, self.parent_pool.fixed_residues):
            individual_protein.residues[ires] = residue
        # put back flipper residues
        for i, residue in zip(self.parent_pool.index_flipper, self.chromosome):
            individual_protein.residues[i] = residue


        individual_protein.prepend_lines = [f"# Fitness = {self.fitness:.2f}; Rank = {self.rank:d}\n"]
        individual_protein.dump(fname)

    def test_lineage(self):  # for debugging
        """
        Test the lineage of this , conf[0] should have parent residue as the original residue while conf[1] should have parent residue as the new residue
        """
        for res in self.chromosome:
            print(f"{res._flag} (=GA)-> {res.conformers[0].parent_residue._flag} (=Empty), {res.conformers[1].parent_residue._flag} (=GA)")


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




