"""
Strip off cofactors if they are not buried in the protein.
"""

from .constants import LOOSE_COFACTORS, UNASSIGEDN_RAD
from ._sas import sas_atoms, radius
import logging

def strip_cofactors(self):      # Here, self is a MCCE object
    """
    Strip off cofactors if they are not buried in the protein.
    """
    sas_cutoff = float(self.prm.SAS_CUTOFF.value)
    loose_cofactors = LOOSE_COFACTORS + self.prm._additional_loose_cofactors # load the default list, will expand this by self.prm._addional_loose_cofactors 


    # collect all atoms in the protein as background atoms
    all_atoms = []
    for res in self.protein.residues:
        all_atoms += res.conformers[0].atoms
        if len(res.conformers) > 1:
            all_atoms += res.conformers[1].atoms

    # check if atom radii are loaded
    for atom in all_atoms:
        if atom.r_vdw < 0.0001:  # if the atom radius is not loaded
            atom.r_vdw = radius.get(atom.element, UNASSIGEDN_RAD)
            logging.debug(f"   atom \"{atom.element}\" has no VDW radius, using default {atom.r_vdw}")



    # loop over all residues and identify loose cofactors from the protein
    # - get resid in a set as potential removers
    # - get their atoms as candidates to calculate SAS
    # - calculate SAS of the cofactor atoms
    # - if SAS is larger than the threshold, remove the residue from the set
    exclude_cofactors = set()
    for res in self.protein.residues:
        if res.resname in LOOSE_COFACTORS:
            res_atoms = res.conformers[0].atoms
            if len(res.conformers) > 1:
                res_atoms += res.conformers[1].atoms
            background = list(set(all_atoms) - set(res_atoms))
            sas = sas_atoms(res_atoms, background)
            sas_exposed = sas_atoms(res_atoms, [])
            sas_percentage = sas / sas_exposed
            #print(f"{res.resname} {res.resid} {sas} {sas_exposed} {sas_percentage}")
            if sas_percentage > sas_cutoff:
                exclude_cofactors.add(res.resid)
    

    # remove the cofactors from the protein by composing a new residues list for protein without exposed cofactors
    new_residues = [res for res in self.protein.residues if res.resid not in exclude_cofactors]
    self.protein.residues = new_residues

    return len(exclude_cofactors)

