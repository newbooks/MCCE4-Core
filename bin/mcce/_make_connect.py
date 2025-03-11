"""
Handles atom connectivity
"""
import logging
from .pdbio import is_H

_BONDDISTANCE_scaling = 0.54  # calibrated by 1akk


# compose ligand detection rules
def compose_ligand_rule(tpl):
    ligand_rules = {}
    for key in tpl.keys():
        if key[0] == "LIGAND_ID":
            res1_name = tpl[key].res1_name  # after renaming
            res2_name = tpl[key].res2_name  # after renaming
            atom1_name = tpl[key].atom1
            atom2_name = tpl[key].atom2
            distance = tpl[key].distance
            tolerance = tpl[key].tolerance
            ligand_rules[(res1_name, res2_name)] = (atom1_name, atom2_name, distance, tolerance)
            ligand_rules[(res2_name, res1_name)] = (atom2_name, atom1_name, distance, tolerance)
    return ligand_rules




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


def match_rule2string(rule, string): # string match that supports wildcard "*"
    return all(r == "*" or r == s for r, s in zip(rule, string))


def make_connect12(self):  # Here, self is a MCCE object
    """
    Make 12 connectivity
    This is a complex task. The connectivity is determined by both the CONNECT records and atom distances.
    The special cases are:
    - The 12 connectivity of side chain atoms iss confined within the same conformer, except to backbone, terminal residues and ligands.
    - The 12 connectivity of backbone atoms is allowed to go to neighboring residues.
    - The 12 connectivity of terminal residues is allowed to go to neighboring residues.
    """
    ligand_rules = compose_ligand_rule(self.tpl)
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
                        # NTR case " CA " to  " C  " and " CB " in the next residue
                        if res.resname in {"NTR", "NTG"} and atom.atomname == " CA ":  # we are looking for " C  " and " CB " in the next residue
                            for conf2 in self.protein.residues[i_res+1].conformers:
                                connected_atom = get_atom_by_name(conf2, " C  ")
                                if connected_atom and connected_atom not in atom.connect12:  # avoid duplicates when checking the next unnamed atom
                                    atom.connect12.append(connected_atom)
                                    logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the next residue")
                                    found = True
                                if not found:
                                    connected_atom = get_atom_by_name(conf2, " CB ")
                                    if connected_atom and connected_atom not in atom.connect12: # avoid duplicates when checking the next unnamed atom
                                        atom.connect12.append(connected_atom)
                                        logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the next residue")
                                        found = True
                        # CTR case: " C  " to  " CA " in the previous residue
                        elif res.resname == "CTR" and atom.atomname == " C  ":
                            conf2 = self.protein.residues[i_res-1].conformers[0]
                            connected_atom = get_atom_by_name(conf2, " CA ")
                            if connected_atom:
                                atom.connect12.append(connected_atom)
                                logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the previous residue")
                                found = True

                        else: # all other ligand cases, unnamed atom to unnamed atom in another residue
                            for res2 in self.protein.residues:
                                if res != res2: # don't look in the same residue
                                    for conf2 in res2.conformers:
                                        for atom2 in conf2.atoms:
                                            key = ("CONNECT", atom2.atomname, conf2.conftype)
                                            if key in self.tpl:
                                                connected_atomnames2 = [x.strip() for x in self.tpl[key].connected] # stripped atom names
                                                if "?" in connected_atomnames2:  # atom2 is an eligible connected atom, check distance
                                                    r = (atom.r_vdw + atom2.r_vdw) * _BONDDISTANCE_scaling  # threshod bond distance
                                                    d = atom.xyz.distance(atom2.xyz)
                                                    if d < r and atom2 not in atom.connect12:
                                                        atom.connect12.append(atom2)
                                                        found = True
                                                        logging.debug(f"Atom {atom.atomname} is connected to {atom2.atomname} in another residue")
                                                    else:  # use ligand rules to detect more connections
                                                        res1_name = res.resname
                                                        res2_name = res2.resname
                                                        atom1_name = atom.atomname
                                                        atom2_name = atom2.atomname
                                                        key = (res1_name, res2_name)
                                                        if key in ligand_rules:
                                                            atom1_name_inrule, atom2_name_inrule, distance, tolerance = ligand_rules[key]
                                                            #print(atom1_name_inrule, atom2_name_inrule, distance, tolerance, d)
                                                        else:
                                                            continue
                                                        # now match the atom names
                                                        print(f"{atom1_name_inrule} =?= {atom1_name} and {atom2_name_inrule} =?= {atom2_name}")
                                                        if match_rule2string(atom1_name_inrule, atom1_name) and match_rule2string(atom2_name_inrule, atom2_name):
                                                            print(f"Matched ligand rule: {key} {ligand_rules[key]}")
                                                            if distance - tolerance < d < distance + tolerance and atom2 not in atom.connect12:
                                                                atom.connect12.append(atom2)
                                                                found = True
                                                                logging.debug(f"Atom {atom.atomname} is connected to {atom2.atomname} by Ligand Rule in another residue")

                                    if found: # one "?" for one ligand
                                        break 

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
                        # 4.1) connects to NTR/NTG, " CB " to previous residue " CA "
                        if not found:
                            if atom.atomname == " CB " and connected_atomname == " CA ":
                                for conf2 in self.protein.residues[i_res-1].conformers[1:]:  # an assumption is NTR is before this residue
                                    connected_atom = get_atom_by_name(conf2, connected_atomname)
                                    if connected_atom:
                                        atom.connect12.append(connected_atom)
                                        logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the previous residue")
                                        found = True
                        # 4.2) connects to NTR, " C  " to previous residue " CA "
                        if not found:
                            if atom.atomname == " C  " and connected_atomname == " CA ":  # an assumption is NTR is before this residue
                                for conf2 in self.protein.residues[i_res-1].conformers[1:]:
                                    connected_atom = get_atom_by_name(conf2, connected_atomname)
                                    if connected_atom:
                                        atom.connect12.append(connected_atom)
                                        logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the previous residue")
                                        found = True
                        # 5) CTR case, " CA " connects " C  " of CTR
                        if not found:
                            if atom.atomname == " CA " and connected_atomname == " C  ":
                                if i_res < len(self.protein.residues):  # didn't reach the end of the protein
                                    for conf2 in self.protein.residues[i_res+1].conformers[1:]:
                                        if conf2.conftype[:3] == "CTR":
                                            connected_atom = get_atom_by_name(conf2, connected_atomname)
                                            if connected_atom:
                                                atom.connect12.append(connected_atom)
                                                logging.debug(f"Atom {atom.atomname} is connected to {connected_atomname} in the next residue")
                                                found = True

                    if not found:
                        logging.debug(f"Atom {connected_atomname} not found in conformer {conf.confid}")


def print_connect12(self):  # Here, self is a MCCE object
    """
    Print 12 connectivity
    """
    for res in self.protein.residues:
        for conf in res.conformers:
            for atom in conf.atoms:
                key = ("CONNECT", atom.atomname, conf.conftype)
                if key in self.tpl:
                    connected_atomnames = self.tpl[key].connected
                else:
                    connected_atomnames = []
                print(f"{atom.residue_id()} {atom.atomname} {connected_atomnames}")
                for connected_atom in atom.connect12:
                    print(f"   {connected_atom.residue_id()} \'{connected_atom.atomname}\'")
                print()