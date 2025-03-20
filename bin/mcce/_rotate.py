"""
Rotamer making module
"""

import logging
from collections import defaultdict
from .pdbio import *

class RotateRule:
    def __init__(self):
        self.bond = ()  # (atom1, atom2)
        self.affected_atoms = []  # list of atoms that will be affected by the rotation

    def __str__(self):
        return f"{self.bond} -> {self.affected_atoms}"
    
    def __expr__(self):
        return self.__str__()

def prepare_rotate_rules(self):  # Hele self is a MCCE object
    """
    Prepare rotamer making rules.
    The combined rules from two sources: ROTATE rules in ftpl and H atom freedom test
    """
    logging.info("   Prepare rotamer making rules ...")
    # self.rotate_rules is a dictionary, in which key is conftype and value is a list of RotateRule objects
    self.rotate_rules = defaultdict(list)

    # Make rotate rules from ftpl
    for key, value in self.tpl.items():
        if key[0] == "ROTATE":
            resname = key[1]
            rotatble_bonds = value.rotatables
            conftypes = self.tpl[("CONFLIST", resname)]
            # search affected atoms of the side chain confs only
            for conftype in conftypes:
                if conftype[-2:] == "BK":
                    continue
                for bond in rotatble_bonds:
                    rule = RotateRule()
                    rule.bond = bond
                    query_key = ("CONNECT", bond[1], conftype)
                    if query_key in self.tpl:
                        rule.affected_atoms.extend([a for a in self.tpl[query_key].connected if a not in bond and ("?" not in a)])
                        # continue to extend the affected atom list by looking for connected atoms of affected atoms
                        for atom in rule.affected_atoms:  # rule.affected_atoms will be extended in the loop
                            query_key = ("CONNECT", atom, conftype)
                            if query_key in self.tpl:
                                rule.affected_atoms.extend([a for a in self.tpl[query_key].connected if (a not in rule.affected_atoms) and (a not in bond) and ("?" not in a)])
                    else:
                        logging.warning(f"   {resname} {conftype} {bond} not found in CONNECT")
                    self.rotate_rules[conftype].append(rule)
    
    # Make rotate rules from H atom freedom test
    for key, value in self.tpl.items():
        if key[0] == "CONNECT" and key[2][-2:] != "BK":  # only side chain confs
            atomname = key[1]
            conftype = key[2]
            if is_H(atomname):  # find H, then find the bonded heavy atom
                for atom in value.connected:  # this should be the only heavy atom connected to the H
                    if not is_H(atom):
                        query_key = ("CONNECT", atom, conftype)
                        orbital = self.tpl[query_key].orbital
                        connected_hvatoms = [a for a in self.tpl[query_key].connected if not is_H(a)]
                        break
                # print(f"   {atomname} {conftype} {atom}: {orbital} {connected_hvatoms}")
                if orbital == "sp3" and len(connected_hvatoms) == 1:  # rotatable oribital and bond has freedom
                    rule = RotateRule()
                    rule.bond = (connected_hvatoms[0], atom)  # bond is defined as (atom1, atom2)
                    query_key = ("CONNECT", atom, conftype)
                    if query_key in self.tpl:
                        rule.affected_atoms.extend([a for a in self.tpl[query_key].connected if a not in rule.bond])  # should all be H atoms
                elif orbital == "sp3" and len(connected_hvatoms) == 0:
                    rule = RotateRule()
                    rule.bond = (None, atom)  # this is a special case to account for free rotate objects like HOH and NH4+
                    query_key = ("CONNECT", atom, conftype)
                    if query_key in self.tpl:
                        rule.affected_atoms.extend([a for a in self.tpl[query_key].connected if a not in rule.bond])  # should all be H atoms
                        # print(f"   Free rotate object {conftype} {atom}: Affected atoms {rule.affected_atoms}")
                self.rotate_rules[conftype].append(rule)

    # Clean up the duplicates in the rotate rules defined by rotatable bonds
    for key, value in self.rotate_rules.items():
        recorded_bonds = []
        clean_rules = []
        for rule in value:
            if rule.bond not in recorded_bonds:
                recorded_bonds.append(rule.bond)
                clean_rules.append(rule)
        self.rotate_rules[key] = clean_rules


def print_rotate_rules(self):
    """
    Print the rotate rules
    """
    for key, value in self.rotate_rules.items():
        print(f"{key}:")
        for rule in value:
            print(rule)