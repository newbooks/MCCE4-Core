"""
Genetic Algorithm for MCCE
"""

import logging
import random
from .pdbio import Residue
from .constants import *
from ._make_connect import get_atom_by_name

class Individual:
    def __init__(self):
        self.chromosome = []  # a list of residues
        self.fitness = None  # fitness score for this individual


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
        self.size = size
        self.index_fixed, self.index_flipper = self.divide_fixed_flipper(mcce)
        self.fixed_residues = [mcce.protein.residues[i] for i in self.index_fixed]
        self.population = []
        # create individuals. pay attention to the performance here
        for i in range(size):
            self.population.append(self.create_individual(mcce))
    
    def divide_fixed_flipper(self, mcce):
        index_fixed = []
        index_flipper = []
        for i, residue in enumerate(mcce.protein.residues):
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

    def create_individual(self, mcce):
        """
        Create an individual for the pool
        """
        individual = Individual()
        for i in self.index_flipper:
            residue = mcce.protein.residues[i]
            # randomly select a conformer from the flipper residues
            selected_conformer = random.choice(residue.conformers).clone()
            selected_conformer.history = selected_conformer.history[:2] + "G" + selected_conformer.history[3:]  # mark as GA selected
            residue = Residue()
            residue.resname = mcce.protein.residues[i].resname
            residue.chain = mcce.protein.residues[i].chain
            residue.sequence = mcce.protein.residues[i].sequence
            residue.insertion = mcce.protein.residues[i].insertion
            residue.conformers = [mcce.protein.residues[i].conformers[0], selected_conformer]  # backbone is a reference while selected_conformer is a copy
            for atom in selected_conformer.atoms:
                key = ("CONNECT", atom.atomname, selected_conformer.conftype)
                if key in mcce.tpl:
                    connected_atomnames = mcce.tpl[key].connected
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
            individual.chromosome.append(residue)
            
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
    