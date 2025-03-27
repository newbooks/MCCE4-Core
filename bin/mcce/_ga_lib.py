"""
Genetic Algorithm for MCCE
"""

import os
import random
from collections import deque
from .pdbio import Residue, Protein
from .constants import *
from ._make_connect import get_atom_by_name

class Individual:
    def __init__(self, pool=None):
        self.parent_pool = pool
        self.chromosome = []  # a list of residues
        self.fitness = 0.0  # fitness score for this individual
        self.rank = 0  # rank for this individual (0 based)

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
        # input parameters
        self.mcce = mcce
        self.size = size
        # own properties
        self.index_fixed, self.index_flipper = self.divide_fixed_flipper()
        self.population = []
        self.pfa_queue = deque(maxlen=GA_PFA_queue)  # keep the last pfa values
        # create individuals. pay attention to the performance here
        self.mcce.reset_connect()  # reset all connect12 and connect13
        self.mcce.assign_qr()  # assign charges and radii
        for i in range(size):
            self.population.append(self.create_individual())
    
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

    def create_individual(self):
        """
        Create an individual for the pool
        """
        individual = Individual(pool=self)
        for i in self.index_flipper:
            residue = self.mcce.protein.residues[i]
            # randomly select a conformer from the flipper residues
            selected_conformer = random.choice(residue.conformers[1:]).clone()
            selected_conformer.history = selected_conformer.history[:2] + "G" + selected_conformer.history[3:]  # mark as GA selected
            new_residue = Residue()
            new_residue.resname = residue.resname
            new_residue.chain = residue.chain
            new_residue.sequence = residue.sequence
            new_residue.insertion = residue.insertion
            new_residue._flag = "GA"
            new_residue.conformers = [residue.conformers[0], selected_conformer]  # backbone is a reference while selected_conformer is a copy
            selected_conformer.parent_residue = new_residue  # note conf[0] still maintain its parent residue to the original residue
            for atom in selected_conformer.atoms:
                key = ("CONNECT", atom.atomname, selected_conformer.conftype)
                if key in self.mcce.tpl:
                    connected_atomnames = self.mcce.tpl[key].connected
                else:
                    connected_atomnames = []
                for connected_atomname in connected_atomnames:
                    connected_atom = get_atom_by_name(selected_conformer, connected_atomname)
                    if connected_atom:
                        atom.connect12.append(connected_atom)
                    else:
                        connected_atom = get_atom_by_name(residue.conformers[0], connected_atomname)
                        if connected_atom:
                            atom.connect12.append(connected_atom)
                # print(atom.atomname, [a.atomname for a in atom.connect12])
            # make connect13
            for atom in selected_conformer.atoms:  # again, we only need this for side chains
                atom.connect13 = []
                for connected_atom in atom.connect12:
                    for connected_atom2 in connected_atom.connect12:
                        if connected_atom2 not in atom.connect13 and connected_atom2 not in atom.connect12 and connected_atom2 != atom:
                            atom.connect13.append(connected_atom2)
                # print(atom.atomname, [a.atomname for a in atom.connect12], [a.atomname for a in atom.connect13])
            individual.chromosome.append(new_residue)     
        return individual

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




