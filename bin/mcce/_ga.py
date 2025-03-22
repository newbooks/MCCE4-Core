"""
Genetic Algorithm for MCCE
"""

import logging
import random
from .pdbio import Residue, Protein
from .constants import *
from ._make_connect import get_atom_by_name
        

class Individual:
    def __init__(self, pool=None):
        self.parent_pool = pool
        self.chromosome = []  # a list of residues
        self.fitness = None  # fitness score for this individual

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


        individual_protein.prepend_lines = ["# This is a mccepdb converted from an individual in GA pool\n"]
        individual_protein.dump(fname)

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



    def __init__(self, mcce, size, ph):
        # input parameters
        self.mcce = mcce
        self.size = size
        self.ph = ph
        # own properties
        self.index_fixed, self.index_flipper = self.divide_fixed_flipper()
        self.fixed_residues = [mcce.protein.residues[i] for i in self.index_fixed]
        self.population = []
        self.pfa = None  # Population Fitness Average for the top 50% individuals
        # create individuals. pay attention to the performance here
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
            new_residue.conformers = [residue.conformers[0], selected_conformer]  # backbone is a reference while selected_conformer is a copy
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
            individual.chromosome.append(new_residue)     
        return individual



def ga_optimize(self):  # Here self is an instance of MCCE
    """
    Genetic Algorithm for MCCE
    """
    # gather constants and parameters for GA
    ga_phs = [float(self.prm.GA_PH1.value), float(self.prm.GA_PH2.value), float(self.prm.GA_PH3.value)]
    ga_pool = int(self.prm.GA_POOL.value)
    ga_maxgen = int(self.prm.GA_MAXGEN.value)
    # the following two variables are numbers of crossovers and mutations in one period (generation)
    ga_crossover = GA_crossover
    ga_mutation = GA_mutation

    logging.info(f"      GA pool size:{ga_pool}")
    logging.info(f"      GA maximum generations:{ga_maxgen}")
    logging.info(f"      GA solution pHs:{ga_phs}")
    for ph in ga_phs:
        logging.info(f"      Prepare GA pool at pH {ph:.2f}, may take a while ...")
        pool = Pool(mcce=self, size=ga_pool, ph=ph)
    
    pool.population[0].to_mccepdb("ga_best.pdb")