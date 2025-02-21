"""
Module Name: pdbio

Description:
This module provides MCCE protein and parameter data structure and input/output functions.
"""

from geom import Vector

class Atom:
    """
    Atom class
    """
    def __init__(self):
        # defined by pdb file
        self.atomname = ""          # atom name
        self.altloc = ""            # alternate location indicator
        self.resname = ""           # residue name
        self.chain = ""             # chain ID
        self.sequence = 0           # residue sequence number
        self.insertion = ""         # insertion code
        self.xyz = Vector()         # coordinates
        # extended attributes
        self.r_boundary = 0.0       # boundary radius
        self.charge = 0.0           # charge
        self.r_vdw = 0.0            # van der Waals radius
        self.e_vdw = 0.0            # van der Waals energy well depth
        self.element = ""           # element name
        self.conn12 = []            # list of atoms that are 1-2 bonded
        self.conn13 = []            # list of atoms that are 1-3 bonded
        self.conn14 = []            # list of atoms that are 1-4 bonded
        self.parent_conf = None     # parent conformer


class Conformer:
    """
    Conformer class
    """
    def __init__(self):
        self.confid = ""            # conformer name, unique ID in the protein. resname+chain+sequence+insertion+confnum
        self.altloc = ""            # alternate location indicator
        self.resname = ""           # residue name
        self.chain = ""             # chain ID
        self.sequence = 0           # residue sequence number
        self.insertion = ""         # insertion code
        self.conftype = ""          # conformer type, as defined in the ftpl file
        self.confnum = 0            # conformer number
        self.atoms = []             # list of atoms in the conformer
        self.parent_residue = None  # parent residue
        self.occ = 0.0              # occupancy
        self.history = ""           # history string
        self.charge = 0.0           # net charge
        self.calculated = False     # flag for calculated conformer


class Residue:
    """
    Residue class
    """
    def __init__(self):
        self.resname = ""           # residue name
        self.chain = ""             # chain ID
        self.sequence = 0           # residue sequence number
        self.insertion = ""         # insertion code
        self.conformers = []        # list of conformers in the residue
        self.resid = ()             # residue ID, (resname, chain, sequence, insertion)


class Protein:
    """
    Protein class
    """
    def __init__(self):
        self.residues = []          # list of residues in the protein


