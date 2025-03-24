"""
Genetic Algorithm for MCCE
"""

import os
import logging
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
        # inituialize the sas for atoms and side chanin conformers

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
    - size: number of individuals in the pool
    - index_fixed: indices of fixed residues
    - index_flipper: indices of flipper residues
    - fixed_residues: list of fixed residues. This is shared by all individuals in the pool
    - flipper_residues = individual: list of flipper residues. This list makes an individual
    - population: list of individuals in the pool
    """
    def __init__(self, mcce, size):
        # input parameters
        self.mcce = mcce
        self.size = size
        # own properties
        self.index_fixed, self.index_flipper = self.divide_fixed_flipper()
        self.fixed_residues = [mcce.protein.residues[i] for i in self.index_fixed]
        self.population = []
        self.pfa_queue = deque(maxlen=GA_PFA_queue)  # keep the last pfa values
        # create individuals. pay attention to the performance here
        self.mcce.reset_connect()  # reset all connect12 and connect13
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
            else:  # more than one side chain
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




def ga_optimize(self, writepdb=False):  # Here self is an instance of MCCE
    """
    Genetic Algorithm for MCCE
    """
    # gather constants and parameters for GA
    ga_pool = int(self.prm.GA_POOL.value)
    ga_maxgen = int(self.prm.GA_MAXGEN.value)
    # the following two variables are numbers of crossovers and mutations in one period (generation)
    ga_crossover = GA_crossover
    ga_mutation = GA_mutation

    # Save the current handlers
    logger = logging.getLogger()
    original_handlers = logger.handlers.copy()

    # Create a file handler
    file_handler = logging.FileHandler(GA_PROGRESS, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s", datefmt='%Y-%m-%d %H:%M:%S'))

    logger.handlers = [file_handler]

    # Start the GA optimization
    logging.info("Genetic Algorithm optimization for MCCE:")
    logging.info(f"   GA pool size:{ga_pool}")
    logging.info(f"   GA maximum generations:{ga_maxgen}")
    logging.info(f"   Prepare GA pool ...")
    pool = Pool(mcce=self, size=ga_pool)
    logging.info(f"      Start to evolve ...")
    pool.evolve()
    logging.info(f"      Done evolution.")
    if writepdb:
        logging.info(f"      Write individuals to {GA_OUTPUT_FOLDER}")
        pool.writepdb()
    logging.info("Genetic Algorithm optimization completed.")

    # pool.population[2].test_lineage()
    # Restore the original handlers
    logger.handlers = original_handlers