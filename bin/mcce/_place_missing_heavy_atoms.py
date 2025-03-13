"""
Place missing *heavy atoms* in the structure. All H atoms are ignored.
"""

import glob
import logging
from .pdbio import *


def load_atom_templates(folder):
    """
    Load atom templates from a folder
    It returns a dictionary of atom templates, with residue name as the first key, atom name as the second key, and atom coordinates as the value
    """
    templates = {}
    for fn in glob.glob(f"{folder}/*.pdb"):
        with open(fn) as f:
            for line in f:
                if line.startswith(("HETATM", "ATOM  ")):
                    atom = Atom()
                    atom.load_pdbline(line)
                    if atom.element != " H":
                        templates.setdefault(atom.resname, {})[atom.atomname] = atom.xyz
    return templates


def place_missing_heavy_atoms(self):    # self is a MCCE object
    # connect12 should already be initialized when entering this function
    n_placed = 0

    self._dist_path = os.path.abspath(os.path.join(__file__, "../../.."))
    ideal_structures = os.path.join(self._dist_path, IDEAL_STRUCTURES)
    templates = load_atom_templates(ideal_structures)

    for res in self.protein.residues:
        for conf in res.conformers:
            completed_missing = set()
            for atom in conf.atoms:
                if atom.element != " H":
                    should_be_connected = set([a for a in self.tpl["CONNECT", atom.atomname, conf.conftype].connected if not is_H(a)])
                    # print(atom.atomname, should_be_connected)
                    actually_connected = set([a.atomname for a in atom.connect12])  # at this point, H atoms already removed
                    # print(atom.atomname, actulayy_connected)
                    missing = should_be_connected - actually_connected - {" ?  "} - completed_missing
                    if missing:
                        n_connected = len(actually_connected)
                        if n_connected == 2:  # 2 known connected atoms, place the missing atom
                            resname = atom.resname
                            



    return n_placed