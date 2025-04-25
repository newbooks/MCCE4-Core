"""
Fast force field calculation for GA.
"""

"""
Fast forcefields for MCCE
Primarily for internal use in Genetic Algorithm optimization.
VDW, Coulomb, and SAS
"""

EPSILON_INNER = 8.0
EPSILON_OUTER = 80.0

import numpy as np

# atom enbedding depth
# Given an indidual in a pool, calculate the enbedding depth of each atom
def atom_embedding_depth(self):  # self is individual in the pool
    """
    Calculate the embedding depth of each atom in the individual.
    The embedding depth is defined as the percentage of 1 A cells that are occupied by atoms.
    """
    # Get all the atoms in the individual
    background_atoms = []
    # fixed atoms from the parent pool
    for i in self.parent_pool.index_fixed:
        for conf in self.parent_pool.mcce.protein.residues[i].conformers:
            background_atoms += conf.atoms
    # atoms from the conf[0] of flipper residues
    for i in self.parent_pool.index_flipper:  # remember, conf[0] keeps lineages to the original residue
        background_atoms += self.parent_pool.mcce.protein.residues[i].conformers[0].atoms
    # get individual (flippable residue) atoms
    individual_atoms = []
    for res in self.chromosome:
        individual_atoms += res.conformers[1].atoms

    # get the dimension of the protein
    coordinates = np.array([[atom.xyz.x, atom.xyz.y, atom.xyz.z] for atom in background_atoms + individual_atoms])
    x_min, y_min, z_min = coordinates.min(axis=0)
    x_max, y_max, z_max = coordinates.max(axis=0)

    # grid the coordinates
    x_bins = np.arange(x_min, x_max + 1, 1)



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







