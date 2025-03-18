"""
Make conformers based on conftype (from CONFLIST)
This module provides the functions to make conformers based on conftype (from CONFLIST).
Prerequsites:
    - MCCE object
    - Heavy atoms are complete
Notes:
    - The hydrogen atoms are not placed at optimized positions
    - The conformers are subjected to further optimization by Genetic Algorithm
"""

import logging
from collections import defaultdict
from .pdbio import *

def propogate_conftypes(self):  # Here self is a MCCE object
    """
    Propogate conftype conformers from CONFLIST
    1. Generate conformers based on conftype
    2. Add hydrogen atoms
    3. Assign charge and radius
    """

    logging.info("   Propogate conftype conformers ...")
    # make a dictionary of conftype and its atom countes (CONNECT records)
    conftype_atom_count = {}
    for key, value in self.tpl.items():
        if key[0] == "CONNECT":
            conftype = key[2]
            if conftype not in conftype_atom_count:
                conftype_atom_count[conftype] = 0
            conftype_atom_count[conftype] += 1
    
    # make conformers based on conftype
    for res in self.protein.residues:
        if ("CONFLIST", res.resname) in self.tpl:
            new_confs = []
            conftype_names = self.tpl[("CONFLIST", res.resname)]
            print(res.resname, conftype_names)
