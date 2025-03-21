"""
Genetic Algorithm for MCCE
"""

import logging
import random
from .constants import *

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
        individual = []
        for i in self.index_flipper:
            residue = mcce.protein.residues[i]
            # randomly select a conformer from the flipper residues
            selected_conformer = random.choice(residue.conformers)
            selected_conformer.history = selected_conformer.history[:2] + "G" + selected_conformer.history[3:]  # mark as GA selected
            individual.append(selected_conformer.clone())
            
        # mix with fixed residues and create connect12 for each atom
        return individual


def ga_optimize(self):  # Here self is an instance of MCCE
    """
    Genetic Algorithm for MCCE
    """
    # gather constants and parameters for GA
    ga_phs = [float(self.prm.GA_PH1.value), float(self.prm.GA_PH2.value), float(self.prm.GA_PH3.value)]
    ga_pool = int(self.prm.GA_POOL.value)
    ga_maxgen = int(self.prm.GA_MAXGEN.value)
    ga_preserve = GA_preserve
    ga_single_cross = GA_single_cross
    ga_double_cross = GA_double_cross
    ga_mutation = GA_mutation
    ga_weedout = GA_weedout
    logging.info(f"      GA pool size:{ga_pool}")
    logging.info(f"      GA maximum generations:{ga_maxgen}")
    logging.info(f"      GA solution pHs:{ga_phs}")
    logging.info(f"      GA preserve population:{ga_preserve}")
    logging.info(f"      GA single_cross rate:{ga_single_cross}")
    logging.info(f"      GA double_cross rate:{ga_double_cross}")
    logging.info(f"      GA mutation rate:{ga_mutation}")
    logging.info(f"      GA weedout rate:{ga_weedout}")
    logging.info(f"      Prepare GA pool, may take a while ...")
    pool = Pool(mcce=self, size=ga_pool)


