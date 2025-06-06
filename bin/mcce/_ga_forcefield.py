"""
Fast force field calculation for GA.
"""

"""
Fast forcefields for MCCE
Primarily for internal use in Genetic Algorithm optimization.
VDW, Coulomb, and SAS
"""
import numpy as np
from .constants import PROBE_RAD

EPSILON_INNER = 8.0
EPSILON_OUTER = 80.0
GRID_SIZE = 1.0   # size of the grid in Angstrom
GRID_DEPTH = 12.0  # depth to look into, this number should be larger than 2 x (the largest atom radius + the probe radius)


class Grid_Pool:     # self is a pool
    """
    A class to represent a grid system that covers a GA pool.
    """
    def __init__(self):
        pass

    def xyz2index(xyz, x_min, y_min, z_min):
        return


    def index2xyz(index, x_min, y_min, z_min):
        return




# atom embedding depth
# Given an indidual in a pool, calculate the enbedding depth of each atom
def atom_embedding_depth(self):  # self is an individual in a pool
    """
    Calculate the embedding depth of each atom in the individual.
    The embedding depth is defined as the percentage of 1 A cells that are occupied by atoms.
    """
    pass

# VDW
def vdw_atom(atom1, atom2):
    pass

def vdw_conformer(conf1, conf2):
    pass

def vdw_protein(protein):
    pass


# ELE - coulomb E = (332 * q1 * q2) / (r) in unit of kcal/mol
def coulomb_atom(atom1, atom2):
    pass

def coulomb_conformer(conf1, conf2):
    pass

def coulomb_protein(protein):
    pass

def ele_atom(atom1, atom2):
    """
    Coulomb interaction between two atoms modified by atom sas.
    - if both atoms are buried, use inner dielectric constant
    - if both atoms are exposed, use outer dielectric constant
    - if one atom is buried and the other is exposed, use linear interpolation
    - atom burial is a modified the product of atom burial and its aprent conformer burial to account for the proximity to surface
    - the final formula is:
        E = (332 * q1 * q2) / (r) * (epsilon_inner * (b1 * b2) + epsilon_outer * (1 - b1 * b2))
    where b1 and b2 are the burial of atom1 and atom2, respectively.
        b1 = buriality of atom1 * buraility of its parent conformer
        b2 = buriality of atom2 * buraility of its parent conformer
    """
    pass


# RXN
def rxn_atom(atom1, atom2):
    pass

def rxn_conformer(conf1, conf2):
    pass

def rxn_protein(protein):
    pass


# SAS - this alone is not an energy term, but used as a modifier for ELE and RXN
def sas_atom(atom):
    pass

def sas_conformer(conf):
    pass







