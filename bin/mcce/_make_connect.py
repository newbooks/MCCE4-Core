"""
Handles atom connectivity
"""
import logging
from .pdbio import is_H

def get_atom_by_name(conf, atomname):
    """
    Get an atom by atomname from a conformer
    """
    for atom in conf.atoms:
        if atom.atomname == atomname:
            return atom
    return None


def reset_connect12(self):  # Here, self is a MCCE object
    """
    Reset 12 connectivity list
    """
    for res in self.protein.residues:
        for conf in res.conformers:
            for atom in conf.atoms:
                atom.connect12 = []


def make_connect12(self):  # Here, self is a MCCE object
    """
    Make 12 connectivity
    This is a complex task. The connectivity is determined by both the CONNECT records and atom distances.
    The special cases are:
    - The 12 connectivity of side chain atoms iss confined within the same conformer, except to backbone, terminal residues and ligands.
    - The 12 connectivity of backbone atoms is allowed to go to neighboring residues.
    - The 12 connectivity of terminal residues is allowed to go to neighboring residues.
    """
    self.reset_connect12()
    for i_res in range(len(self.protein.residues)): # Loop over residues and we need index number to look up neighboring residues
        res = self.protein.residues[i_res]
        for conf in res.conformers:
            for atom in conf.atoms:
                key = ("CONNECT", atom.atomname, conf.conftype)
                if key in self.tpl:
                    connected_atomnames = self.tpl[key].connected
                    logging.debug(f"Atom {atom.atomname} is connected to {connected_atomnames}")
                else:
                    connected_atomnames = []
                for connected_atomname in connected_atomnames:
                    logging.debug(f"Looking for connected atom {connected_atomname}")
                    found = False   # This is a flag to indicate if the connected atom is found
                    if "?" in connected_atomname:  # "?" mean unnamed atom, this is a connected atom outside the residue
                        pass
                    
                    else:  # named atom, within the same conformer, or from side chain to the backbone, or from the backbbone to all side chains
                        # 1) within the same conformer
                        connected_atom = get_atom_by_name(conf, connected_atomname)
                        if connected_atom: # found within the same conformer, all good!
                            atom.connect12.append(connected_atom)
                            logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the same conformer")
                            found = True
                        if not found: # not found within the same conformer, 
                            # 2) from backbone atome, check all side chains
                            if conf.conftype[-2:] == "BK":  # if conf is backbone, look for the atom in all side chains
                                for conf2 in res.conformers:
                                    if conf2.conftype[-2:] != "BK":
                                        connected_atom = get_atom_by_name(conf2, connected_atomname)
                                        if connected_atom:
                                            atom.connect12.append(connected_atom)
                                            logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in side chain")
                                            found = True
                            # 3) from side chain to the backbone
                            else:  # if conf is side chain, look for the atom in the backbone
                                connected_atom = get_atom_by_name(res.conformers[0], connected_atomname)
                                if connected_atom:
                                    atom.connect12.append(connected_atom)
                                    logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the backbone")
                                    found = True
                        # 4) connects to NTR/NTG, " CB " to previous residue " CA "
                        if not found:
                            if atom.name == " CB " and connected_atomname == " CA ":
                                for conf2 in self.protein.residues[i_res-1].conformers[1:]:
                                    connected_atom = get_atom_by_name(conf2, connected_atomname)
                                    if connected_atom:
                                        atom.connect12.append(connected_atom)
                                        logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the previous residue")
                                        found = True

                    if not found:
                        logging.debug(f"Atom {connected_atomname} not found in conformer {conf.confid}")
