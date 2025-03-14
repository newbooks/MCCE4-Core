"""
Place missing *heavy atoms* in the structure. All H atoms are ignored.
"""

import sys
import glob
import logging
from .pdbio import *
from .geom import *
from .constants import *

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
        if len(res.conformers) < 2:  # only backbone
            continue
        for conf in res.conformers[1:]: # skip the first conformer, which is the backbone
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
                        if n_connected >= 2:  # at least 2 known connected atoms, place the missing atom
                            resname = atom.resname
                            c_atom_p = atom.xyz         # center atom coordinates
                            c_atom_t = templates[resname][atom.atomname]  # center atom of the template coordinates

                            # will get the two connected atoms
                            check_atom = atom.connect12[0]  # pick one of the connected atoms
                            if check_atom.atomname not in templates[resname]:
                                logging.error(f"   Connected atom \'{check_atom.atomname}\' to \'{atom.atomname}\' not found in template for residue {resname}")
                                sys.exit(1)
                            a_atom_p = check_atom.xyz
                            a_atom_t = templates[resname][check_atom.atomname]

                            check_atom = atom.connect12[1]  # pick the other connected atom
                            if check_atom.atomname not in templates[resname]:
                                logging.error(f"   Connected atom \'{check_atom.atomname}\' to \'{atom.atomname}\' not found in template for residue {resname}")
                                sys.exit(1)
                            b_atom_p = check_atom.xyz
                            b_atom_t = templates[resname][check_atom.atomname]

                            # now we have all the coordinates, we can place the missing atom
                            op = geom_3v_onto_3v(c_atom_t, a_atom_t, b_atom_t, c_atom_p, a_atom_p, b_atom_p)
                            for atom_name in list(missing):
                                new_atom = Atom()
                                new_atom.atomname = atom_name
                                new_atom.xyz = op.apply_to_vector(templates[resname][atom_name])  # apply the transformation to the template coordinates
                                conf.atoms.append(new_atom)
                                completed_missing.add(atom_name)

    return n_placed