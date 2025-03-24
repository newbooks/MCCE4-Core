"""
Fast forcefields for MCCE
Primarily for internal use in Genetic Algorithm optimization.
VDW, Coulomb, and SAS
"""

EPSILON_INNER = 8.0
EPSILON_OUTER = 80.0

import numpy as np

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







