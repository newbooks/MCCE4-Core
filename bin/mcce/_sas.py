"""
This module contains functions for calculating solvent accessible surface area (SAS) of a protein.
- sas_atom(atom, background): a function that calculates SAS of an atom in context of background.
- sas_residue(residue, background): a function that calculates SAS of a residue in context of background.

SAS is calculated using van der Waals radii of atoms and a rolling ball algorithm.
"""

import math
import os
from collections import defaultdict
import logging

from .geom import *
from .constants import *
from .pdbio import *

# default radius of atoms, will be overwritten by non-0 r_vdw in atom
radius = {" H": 1.2,
          " C": 1.7,
          " N": 1.55,
          " O": 1.52,
          " F": 1.47,
          " P": 1.8,
          " S": 1.8,
          "CL": 1.75,
          "CU": 1.4,
          " B": 1.92,
          "AL": 1.84,
          "NA": 2.27,
          "MG": 1.73,
          "SI": 2.1,
          "CA": 2.31,
          " K": 2.75,
          "FE": 1.63,
          "ZN": 1.39,
          "BR": 1.85,
          " ?": UNASSIGEDN_RAD
}

probe_global = PROBE_RAD  # global variable for probe radius, it gets default value from constants.py, can be updated by other subroutines


# preset points on a sphere for SAS calculation
n_points = 89   # sampling points on a sphere, Fibonacci numbers: 34,55,89,144,233...

def fibonacci_sphere(n):
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
    return [Vector((math.cos(phi * i) * math.sqrt(1 - (1 - (i / float(n - 1)) * 2) ** 2),
                      1 - (i / float(n - 1)) * 2,
                      math.sin(phi * i) * math.sqrt(1 - (1 - (i / float(n - 1)) * 2) ** 2)))
              for i in range(n)]

points_preset = fibonacci_sphere(n_points)


def sas_atom(atom, background):
    """
    Compute atom sas in the context of background atoms 
    """

    # get atoms in the background that can have impact on the SAS of the atom
    relavent_background = []
    for b in background:
        max_dist = atom.r_vdw + b.r_vdw + 2 * probe_global  # maximum distance between atom and b that can have impact on SAS of atom
        if abs(b.xyz.x - atom.xyz.x) < max_dist and abs(b.xyz.y - atom.xyz.y) < max_dist and abs(b.xyz.z - atom.xyz.z) < max_dist:
            relavent_background.append(b)
    #print(len(background), len(relavent_background))

    # generate points on a sphere around the atom
    points_on_sphere = [point * (atom.r_vdw + probe_global) + atom.xyz for point in points_preset]
    n_buried = 0
    for point in points_on_sphere:
        for b in relavent_background:
            v = point - atom.xyz
            dd = v.x * v.x + v.y * v.y + v.z * v.z
            if dd < (b.r_vdw + probe_global)**2:
                n_buried += 1
                break
            
    return n_buried / n_points * 4 * math.pi * (atom.r_vdw + probe_global) ** 2

def sas_residue(residue, background):
    """
    Compute residue sas in the context of background atoms
    """

def sas_atoms(atoms, background):
    """
    Compute SAS for a list of atoms
    It returns the total solvent accessible surface area of the atoms, and stores individual SAS in the atom object attribute sas.
    """
    # find sas of each atom
    total_sas = 0
    for atom in atoms:
        other_atoms = list(set(atoms) - {atom})
        atom.sas = sas_atom(atom, background + other_atoms)
        total_sas += atom.sas
    
    return total_sas


def sas_pdb(pdb, probe):
    """
    Compute SAS for a pdb file
    """
    # update global probe radius with the input probe radius
    global probe_global
    probe_global = probe

    # assign radius to atoms
    for atom in pdb.atoms:
        if atom.r_vdw < 0.0001:
            atom.r_vdw = radius.get(atom.element, UNASSIGEDN_RAD)
            logging.debug(f"Warning: atom \"{atom.element}\" has no VDW radius, using default {atom.r_vdw}")

    # group atoms by residue
    residues = defaultdict(list)
    for atom in pdb.atoms:
        residues[atom.residue_id()].append(atom)

    # calculate SAS for each residue
    all_atoms = set(pdb.atoms)
    for res, atoms in residues.items():
        background = list(all_atoms - set(atoms))
        sas_atoms(atoms, background)