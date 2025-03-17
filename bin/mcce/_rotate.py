"""
Rotamer making module
"""

import logging
from collections import defaultdict

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
                        rule.affected_atoms.extend([a for a in self.tpl[query_key].connected if a not in bond])
                        # continue to extend the affected atom list by looking for connected atoms of affected atoms
                        for atom in rule.affected_atoms:  # rule.affected_atoms will be extended in the loop
                            query_key = ("CONNECT", atom, conftype)
                            if query_key in self.tpl:
                                rule.affected_atoms.extend([a for a in self.tpl[query_key].connected if a not in rule.affected_atoms])
                    else:
                        logging.warning(f"   {resname} {conftype} {bond} not found in CONNECT")
                    self.rotate_rules[conftype].append(rule)
    
    for key, value in self.rotate_rules.items():
        print(f"{key}:")
        for rule in value:
            print(rule)