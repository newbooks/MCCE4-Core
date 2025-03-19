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
    #   r3 r4 r5             y
    #     \|/               /
    #      r0           x--0
    #      |               |
    #      r1              z (v10)
    #      \                \  
    #       r2               r2
    # r0 is sp3 type.
    # r0, r1, r2's coordinates are known.
    # r3, r4, r5's coordinates are to be determined.
    #
    # Solution:
    # 1. Calculate uz, this is the direction of r1 from r0
    # 2. Normalize uz, this is the unit vector of z axis
    # 3a. If r2 is not defined, pick an arbitrary orthogonal vector to uz as uy
    # 3b. Calculate v012, this is the direction of plane r0, r1, r2
    #     Calculate uy = v012/|v012|, this is the unit vector of v012
    # 4. Calculate ux = uy x uz, this is the unit vector of uy cross uz
    # 5. ux, uy, uz is the new coordinate system, with uz (r01's unit vector) as z axis and other two as x and y axes
    # 6. Calculate r3, r4, r5 in the new coordinate system
    # 7. return r3, r4, r5 in the original coordinate system
    uz = r1 - r0
    uz.normalize()
    
    if r2 is None:
        uy = uz.orthogonal()
    else:
        u12 = r2 - r1
        uy = uz.cross(u12)
    uy.normalize()

    ux = uy.cross(uz)
    ux.normalize()

    # now we have the new coordinate system
    cos_theta = math.cos(math.radians(H_BOND_ANGLE_SP3))
    sin_theta = math.sin(math.radians(H_BOND_ANGLE_SP3))

    u3 = uz * cos_theta + ux * sin_theta
    u3.normalize()
    r3 = r0 + u3 * H_BOND_LENGTH

    x_component = ux * (sin_theta * cos_theta)
    y_component = uy * (sin_theta * sin_theta)
    z_component = uz * cos_theta
    r4 = r0 + (x_component + y_component + z_component) * H_BOND_LENGTH

    y_component = uy * (- sin_theta * cos_theta)
    r5 = r0 + (x_component + y_component + z_component) * H_BOND_LENGTH
    
    return r3, r4, r5



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
                                if len(missing_H) == 2:
                                    new_atom1 = create_connected_atom(atom, missing_H.pop(0), r3)
                                    new_atom2 = create_connected_atom(atom, missing_H.pop(0), r4)
                                    conf.atoms.append(new_atom1)
                                    conf.atoms.append(new_atom2)
                                elif len(missing_H) == 1:
                                    new_atom = create_connected_atom(atom, missing_H.pop(0), r3)
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