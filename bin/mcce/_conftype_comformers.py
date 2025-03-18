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
from .constants import *
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
    for key, _ in self.tpl.items():
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
            nondummy_conftype_names = set([c for c in conftype_names if c in conftype_atom_count and c[-2:] != "BK"])
            # print(res.resname, conftype_names, nondummy_conftype_names)

            # propogate from the existing non-BK conformers to the new conformers
            if len(res.conformers) > 1:
                for parent_conf in res.conformers[1:]:
                    other_conftypes = nondummy_conftype_names -set([parent_conf.conftype])
                    for conftype in other_conftypes:
                        new_conf = parent_conf.clone()
                        new_conf.conftype = conftype
                        new_conf.history = conftype[-2:] + parent_conf.history[2:]
                        new_confs.append(new_conf)
            
                # merge the new conformers to existing conformers
                all_sidechain_confs = res.conformers[1:] + new_confs
                ordered_confs = []
                for conftype in conftype_names:
                    ordered_confs.extend([c for c in all_sidechain_confs if c.conftype == conftype])  # sort by conftype in CONFLIST
                res.conformers = [res.conformers[0]] + ordered_confs
        else:
            logging.warning(f"   {res.resname} has no CONFLIST record")

    # add hydrogen atoms
    self.make_connect12()  # make connect12 records again as new conformers are added

    for res in self.protein.residues:
        for conf in res.conformers:
            for atom in conf.atoms:
                if atom.element == " H":
                    continue
                should_have_H = [a for a in self.tpl[("CONNECT", atom.atomname, conf.conftype)].connected if is_H(a)]
                actually_have_H = [a.atomname for a in atom.connect12 if is_H(a.atomname)]
                missing_H = [a for a in should_have_H if a not in actually_have_H]
                if missing_H:
                    orbital = self.tpl[("CONNECT", atom.atomname, conf.conftype)].orbital
                    connected_heavy_atoms = [a for a in atom.connect12 if not is_H(a.atomname)]
                    # backbone atom may have multiple side chain heavy atoms connected in old implementation, so we filter it out
                    connected_heavy_atoms_byname = set()
                    uniq_connected_heavy_atoms = []
                    for a in connected_heavy_atoms:
                        if a.atomname not in connected_heavy_atoms_byname:
                            connected_heavy_atoms_byname.add(a.atomname)
                            uniq_connected_heavy_atoms.append(a)
                    connected_heavy_atoms = uniq_connected_heavy_atoms

                    # decide how to place H atom
                    if orbital.lower() == "ion":
                        continue
                    elif orbital.lower() == "sp3":
                        n_known = len(connected_heavy_atoms)
                        if n_known == 3:
                            #    r4
                            #     |
                            #    r0
                            #    /|\
                            # r1 r2 r3
                            #
                            # r0 is sp3 type.
                            # r0, r1, r2, r3's coordinates are known.
                            # r4's coordinate is to be determined.
                            if len(missing_H) == 1:
                                r0 = atom
                                r1 = connected_heavy_atoms[0]
                                r2 = connected_heavy_atoms[1]
                                r3 = connected_heavy_atoms[2]
                                r4 = Atom()
                                r4.atomname = missing_H[0]
                                r4.resname = res.resname
                                r4.resid = res.resid
                                r4.chainid = res.chainid
                                r4.xyz = get_H_xyz(r0, r1, r2, r3)
                                r4.history = "H"
                                res.add_atom(r4)
