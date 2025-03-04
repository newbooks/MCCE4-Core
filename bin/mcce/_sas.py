"""
This module contains functions for calculating solvent accessible surface area (SAS) of a protein.
- sas_atom(atom, background): a function that calculates SAS of an atom in context of background.
- sas_residue(residue, background): a function that calculates SAS of a residue in context of background.

SAS is calculated using van der Waals radii of atoms and a rolling ball algorithm.
"""

import math
import os
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

box_expand = max(radius.values()) + PROBE_RAD   # expand the box by the largest atom radius and probe radius


def fibonacci_sphere(n):
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
    return [Vector((math.cos(phi * i) * math.sqrt(1 - (1 - (i / float(n - 1)) * 2) ** 2),
                      1 - (i / float(n - 1)) * 2,
                      math.sin(phi * i) * math.sqrt(1 - (1 - (i / float(n - 1)) * 2) ** 2)))
              for i in range(n)]


def sas_atom(atom, background):
    """
    Compute atom sas in the context of background atoms 
    """

def sas_residue(residue, background):
    """
    Compute residue sas in the context of background atoms
    """

def sas_pdb(pdb_file, probe):
    """
    Compute SAS for a pdb file
    """
    if os.path.exists(pdb_file):
        pdb = Pdb(pdb_file)
    else:
        print(f"File {pdb_file} does not exist.")
        return
    