"""
This module contains functions for calculating solvent accessible surface area (SAS) of a protein.
- sas_atom(atom, background): a function that calculates SAS of an atom in context of background.
- sas_residue(residue, background): a function that calculates SAS of a residue in context of background.

SAS is calculated using van der Waals radii of atoms and a rolling ball algorithm.
"""

def sas_atom(atom, background):
    """
    Compute atom sas in the context of background atoms 
    """