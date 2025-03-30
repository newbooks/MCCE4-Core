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
import math
from collections import defaultdict
from .constants import *
from .pdbio import *

atom_mass = {" H" : 1,
             " C" : 12,
                " N" : 14,
                " O" : 16,
                " S" : 32,
                " P" : 31,
                " F" : 19,
                "CL": 35.5,
                "BR": 80,
                "I" : 127,
                "ZN": 65,
                "CA": 40,
                "MG": 24,
                "MN": 55,
                "FE": 56,
                "CU": 63,
                "CO": 59,
                "NI": 59,
                " B": 10,
                "SI": 28,
                "LI": 7,
                "BE": 9,
                "NA": 23,
                "AL": 27,
                " K": 39}


def create_connected_atom(atom, atomname, xyz):
    """
    Create a new atom based on an existing atom
    """
    new_atom = atom.clone()
    new_atom.atomname = atomname # new atom name
    new_atom.xyz = xyz  # new coordinates
    if len(atomname.strip()) == 4 and atomname[0] == "H": # new element
        new_atom.element = " H"
    else:
        new_atom.element = atomname[:2]
    new_atom.connect12.append(atom)  # add the original atom to connect12
    atom.connect12.append(new_atom)  # add the new atom to the original atom's connect12
    return new_atom


def sp3_2knowns(r0, r1, r2):
    """
    Place r3 and r4 connected to r0 when r1 and r2 are known
    """
    #   r3 r4
    #     \/
    #      r0 - r2
    #      |
    #      r1
    #
    # r0 is sp3 type.
    # r0, r1, r2's coordinates are known.
    # r3, r4's coordinates are to be determined.
    #
    # Solution: r3, r4 are placed at the opposite direction of the bisector of r1 and r2
    u1 = r1 - r0
    u2 = r2 - r0
    u1.normalize()
    u2.normalize()
    bisector = (u1 + u2)
    bisector.normalize()
    opposite_bisector = bisector * (-1)  # opposite direction
    norm102 = u1.cross(u2)
    norm102.normalize()
    half_angle = math.radians(H_BOND_ANGLE_SP3) / 2
    cosine = math.cos(half_angle)
    sine = math.sin(half_angle)
    r3 = r0 + (opposite_bisector * cosine + norm102 * sine) * H_BOND_LENGTH
    r4 = r0 + (opposite_bisector * cosine - norm102 * sine) * H_BOND_LENGTH
    return r3, r4

def sp3_1known(r0, r1, r2):
    """
    Place 3 H (r3, r4, and r5) atoms to r0
    """
    #   r3 r4 r5          
    #     \|/             
    #      r0              0 - x
    #      |             / |
    #      r1            y  z (v01)
    #      \                \  
    #       r2               r2
    # r0 is sp3 type.
    # r0, r1, r2's coordinates are known.
    # r3, r4, r5's coordinates are to be determined.
    uz = r1 - r0
    uz.normalize()  # pointing down
    
    if r2 is None:
        u2 = uz.orthogonal()
    else:
        u2 = r2 - r1

    uy = uz.cross(u2)  
    uy.normalize()    # pointing towards us

    ux = uy.cross(uz)
    ux.normalize()  # pointing to the right


    # now we have the new coordinate system
    cos_theta = math.cos(math.radians(H_BOND_ANGLE_SP3))
    sin_theta = math.sin(math.radians(H_BOND_ANGLE_SP3))

    u3 = uz * cos_theta - ux * sin_theta
    u3.normalize()
    r3 = r0 + u3 * H_BOND_LENGTH

    x_component = ux * (-sin_theta * cos_theta)
    y_component = uy * (sin_theta * sin_theta)
    z_component = uz * cos_theta
    r4 = r0 + (x_component + y_component + z_component) * H_BOND_LENGTH

    y_component = uy * (-sin_theta * sin_theta)
    r5 = r0 + (x_component + y_component + z_component) * H_BOND_LENGTH
    
    return r3, r4, r5


def sp3_0known(r0):
    """
    place 4 H (r1, r2, r3, r4) atoms to r0
    """
    #   r1  r2
    #     \ /
    #      r0--r3
    #      |
    #     r4
    #
    # r0 is sp3 type.
    # r0's coordinates are known.
    # r1, r2, r3 and r4's coordinates are to be determined.
    #
    # Solution: r1, r2, r3, r4 are placed at the vertices of a regular tetrahedron with r0 as the center (0,0,0)
    u1 = Vector((1, 1, 1))
    u2 = Vector((-1, -1, 1))
    u3 = Vector((-1, 1, -1))
    u4 = Vector((1, -1, -1))
    u1.normalize()
    u2.normalize()
    u3.normalize()
    u4.normalize()
    r1 = r0 + u1 * H_BOND_LENGTH
    r2 = r0 + u2 * H_BOND_LENGTH
    r3 = r0 + u3 * H_BOND_LENGTH
    r4 = r0 + u4 * H_BOND_LENGTH
    return r1, r2, r3, r4


def sp2_2known(r0, r1, r2):
    """
    Place r3 connected to r0 when r1 and r2 are known
    """
    #   r3
    #     \
    #      r0 - r2
    #      |
    #      r1
    #
    # r0 is sp2 type.
    # r0, r1, r2's coordinates are known.
    # r3's coordinates are to be determined.
    #
    # Solution: r3 is placed at the bisector of r1 and r2
    u1 = r1 - r0
    u2 = r2 - r0
    u1.normalize()
    u2.normalize()
    bisector = u1 + u2
    bisector.normalize()
    r3 = r0 - bisector * H_BOND_LENGTH
    return r3


def sp2_1known(r0, r1, r2):
    """
    Place 2 H (r3 and r4) atoms to r0
    """
    #   r3 r4
    #     \/
    #      r0
    #      |
    #      r1
    #      /
    #     r2
    #
    # r0 is sp2 type.
    # r0, r1, r2's coordinates are known.
    # r3, r4's coordinates are to be determined.
    #
    # Solution: r3, r4 are placed on the same plane defined by r0, r1, and r2, with bond angle of 120 degrees
    uy = r1 - r0
    uy.normalize()
    if r2 is None:
        v12 = uy.orthogonal()
    else:
        v12 = r2 - r1
    uz = uy.cross(v12)
    uz.normalize()
    ux = uy.cross(uz)
    ux.normalize()

    cos_theta = math.cos(math.radians(H_BOND_ANGLE_SP2 - 90))
    sin_theta = math.sin(math.radians(H_BOND_ANGLE_SP2 - 90))
    u4 = ux * cos_theta - uy * sin_theta
    u4.normalize()
    r4 = r0 + u4 * H_BOND_LENGTH
    u3 = ux * (-cos_theta) - uy * sin_theta
    u3.normalize()
    r3 = r0 + u3 * H_BOND_LENGTH
    return r3, r4


def sp2_0known(r0):
    """
    Place 3 H (r1, r2, r3) atoms to r0
    """
    #   r1  r2
    #     \ /
    #      r0
    #      |
    #      r3
    #
    # r0 is sp2 type.
    # r0's coordinates are known.
    # r1, r2, r3's coordinates are to be determined.
    #
    # Solution: r1, r2, r3 are placed at the vertices of an equilateral triangle with r0 as the center (0,0,0)
    angle = math.radians(120)
    u1 = Vector((math.cos(0), math.sin(0), 0))
    u2 = Vector((math.cos(angle), math.sin(angle), 0))
    u3 = Vector((math.cos(-angle), math.sin(-angle), 0))
    r1 = r0 + u1 * H_BOND_LENGTH
    r2 = r0 + u2 * H_BOND_LENGTH
    r3 = r0 + u3 * H_BOND_LENGTH
    return r1, r2, r3


def propogate_conftypes(self):  # Here self is a MCCE object
    """
    Propogate conftype conformers from CONFLIST
    1. Generate conformers based on conftype
    2. Add hydrogen atoms
    3. Assign charge and radius
    """

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
                                r0 = atom.xyz
                                r1 = connected_heavy_atoms[0].xyz
                                r2 = connected_heavy_atoms[1].xyz
                                r3 = connected_heavy_atoms[2].xyz
                                u1 = r1 - r0
                                u2 = r2 - r0
                                u3 = r3 - r0
                                u1.normalize()
                                u2.normalize()
                                u3.normalize()
                                u_sum = u1 + u2 + u3
                                u_sum.normalize()
                                r4 = r0 - u_sum * H_BOND_LENGTH
                                new_atom = create_connected_atom(atom, missing_H.pop(0), r4)
                                conf.atoms.append(new_atom)
                            else:
                                logging.warning(f"   {atom.atomname} {conf.conftype} has {n_known} connected heavy atoms, but {len(missing_H)} missing H atoms. Cannot place H.")
                        elif n_known == 2:
                            # r3 r4
                            #   \/
                            #   r0
                            #   /|
                            # r1 r2
                            #
                            # r0 is sp3 type.
                            # r0, r1, r2's coordinates are known.
                            # r3, r4's coordinates are to be determined.
                            if len(missing_H) <= 2:
                                r0 = atom.xyz
                                r1 = connected_heavy_atoms[0].xyz
                                r2 = connected_heavy_atoms[1].xyz
                                r3, r4 = sp3_2knowns(r0, r1, r2)
                                for r in [r3, r4][:len(missing_H)]:
                                    new_atom = create_connected_atom(atom, missing_H.pop(0), r)
                                    conf.atoms.append(new_atom)
                            else:
                                logging.warning(f"   {atom.atomname} {conf.conftype} has {n_known} connected heavy atoms, but {len(missing_H)} missing H atoms. Cannot place H.")
                        elif n_known == 1:
                            #   r3 r4 r5
                            #     \|/
                            #      r0
                            #      |
                            #      r1
                            #      \
                            #       r2
                            # r0 is sp3 type.
                            # r0, r1, r2's coordinates are known.
                            # r3, r4, r5's coordinates are to be determined.
                            torsion_minimum = False
                            if len(missing_H) <= 3:
                                r0 = atom.xyz
                                r1 = connected_heavy_atoms[0].xyz
                                # pick r2, the heaviest atom connected to r1
                                extended_atoms = list(set(connected_heavy_atoms[0].connect12) - set([atom]))
                                if len(extended_atoms) == 1:  # only one atom connected to r1
                                    r2 = extended_atoms[0].xyz
                                elif len(extended_atoms) == 0:
                                    r2 = None
                                else:
                                    # init atom mass
                                    for a in extended_atoms:
                                        if a.element:
                                            if a.element in atom_mass:
                                                a.mass = atom_mass[a.element]
                                            else:
                                                a.mass = 100  # unknown atom, assume to be a heavy atom
                                        else:
                                            logging.error(f"   Atom {a.atomname} in {a.parent_conf.parent_residue.resname} has no element defined")
                                            sys.exit(1)
                                    # sort by mass
                                    extended_atoms.sort(key=lambda x: x.mass, reverse=True)
                                    r2 = extended_atoms[0].xyz
                                torsion_minimum = True
                                r3, r4, r5 = sp3_1known(r0, r1, r2)
                                for r in [r3, r4, r5][:len(missing_H)]:
                                    new_atom = create_connected_atom(atom, missing_H.pop(0), r)
                                    conf.atoms.append(new_atom)

                                if torsion_minimum:
                                    conf.history = conf.history[:6] + "M" + conf.history[7:]  # mark as torsion minimum
                                    torsion_minimum = False

                            else:
                                logging.warning(f"   {atom.atomname} {conf.conftype} has {n_known} connected heavy atoms, but {len(missing_H)} missing H atoms. Cannot place H.")
                            
                        elif n_known == 0:
                            # this applies to heavy atoms that do not have other heavy atoms connected, like H2O and NH4+
                            r0 = atom.xyz
                            r1, r2, r3, r4 = sp3_0known(r0)
                            for r in [r1, r2, r3, r4][:len(missing_H)]:
                                new_atom = create_connected_atom(atom, missing_H.pop(0), r)
                                conf.atoms.append(new_atom)                                
                    
                    elif orbital.lower() == "sp2":            
                        n_known = len(connected_heavy_atoms)
                        if n_known == 2: # place the 3rd H atom at a defined place
                            # r3
                            #  \
                            #   r0 - r2
                            #   |
                            #   r1
                            #
                            # r0 is sp2 type.
                            # r0, r1, r2's coordinates are known.
                            # r3's coordinates are to be determined.
                            if len(missing_H) == 1:
                                r0 = atom.xyz
                                r1 = connected_heavy_atoms[0].xyz
                                r2 = connected_heavy_atoms[1].xyz
                                r3 = sp2_2known(r0, r1, r2)
                                new_atom = create_connected_atom(atom, missing_H.pop(0), r3)
                                conf.atoms.append(new_atom)
                        elif n_known == 1:  # place the 2nd and 3rd H atoms at defined places
                            # r3 r4
                            #   \/
                            #   r0
                            #   |
                            #   r1
                            #   /
                            #  r2
                            #
                            # r0 is sp2 type.
                            # r0, r1, r2's coordinates are known.
                            # r3, r4's coordinates are to be determined.
                            r0 = atom.xyz
                            r1 = connected_heavy_atoms[0].xyz
                            torsion_minimum = False
                            # pick r2, the heaviest atom connected to r1
                            extended_atoms = list(set(connected_heavy_atoms[0].connect12) - set([atom]))
                            if len(extended_atoms) == 1:  # only one atom connected to r1
                                r2 = extended_atoms[0].xyz
                                torsion_minimum = True
                            elif len(extended_atoms) == 0:
                                r2 = None
                            else:
                                # init atom mass
                                for a in extended_atoms:
                                    if a.element:
                                        if a.element in atom_mass:
                                            a.mass = atom_mass[a.element]
                                        else:
                                            a.mass = 100
                                    else:
                                        logging.error(f"   Atom {a.atomname} in {a.parent_conf.parent_residue.resname} has no element defined")
                                        sys.exit(1)
                                # sort by mass
                                extended_atoms.sort(key=lambda x: x.mass, reverse=True)
                                r2 = extended_atoms[0].xyz
                                torsion_minimum = True
                            r3, r4 = sp2_1known(r0, r1, r2)
                            for r in [r3, r4][:len(missing_H)]:
                                new_atom = create_connected_atom(atom, missing_H.pop(0), r)
                                conf.atoms.append(new_atom)
                            if torsion_minimum:
                                conf.history = conf.history[:6] + "M" + conf.history[7:]
                                torsion_minimum = False
                        elif n_known == 0:
                            r0 = atom.xyz
                            r1, r2, r3 = sp2_0known(r0)
                            for r in [r1, r2, r3][:len(missing_H)]:
                                new_atom = create_connected_atom(atom, missing_H.pop(0), r)
                                conf.atoms.append(new_atom)
                    elif orbital.lower() == "sp":
                        n_known = len(connected_heavy_atoms)
                        if n_known == 1:
                            # r2
                            #  \
                            #   r0
                            #    \
                            #     r1
                            r0 = atom.xyz
                            r1 = connected_heavy_atoms[0].xyz
                            r2 = (r0 - r1).normalize() * H_BOND_LENGTH + r0
                            if len(missing_H) == 1:
                                new_atom = create_connected_atom(atom, missing_H.pop(0), r2)
                                conf.atoms.append(new_atom)
                    else:
                        logging.error(f"   Unknown orbital type {orbital} for {atom.atomname} {conf.conftype}")
                        sys.exit(1)

def propogate_swap(self): # Here self is a MCCE object
    """
    Propogate conformers by ROT_SWAP rules
    """
    for res in self.protein.residues:
        query_key = ("ROT_SWAP", res.resname)
        if query_key in self.tpl:
            swap_rule = self.tpl[query_key]
            swap_from = [a[0] for a in swap_rule.swapables]
            swap_to = [a[1] for a in swap_rule.swapables]
            for conf in res.conformers[1:]:
                new_conf = conf.clone()
                swap_from_atoms = [a for a in new_conf.atoms if a.atomname in swap_from]
                swap_to_atoms = [a for a in new_conf.atoms if a.atomname in swap_to]
                if len(swap_from_atoms) != len(swap_from) or len(swap_to_atoms) != len(swap_to):
                    logging.warning(f"   {res.resname} {conf.conftype} has {len(swap_from_atoms)} swap from atoms and {len(swap_to_atoms)} swap to atoms")
                    continue 
                new_conf.history = conf.history[:2] + "W" + conf.history[3:]
                swapped = False
                for i in range(len(swap_from)):
                    swap_from_atoms[i].xyz, swap_to_atoms[i].xyz = swap_to_atoms[i].xyz, swap_from_atoms[i].xyz
                    swapped = True
            if swapped:
                res.conformers.append(new_conf)