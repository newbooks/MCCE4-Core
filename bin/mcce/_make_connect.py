"""
Handles atom connectivity
"""

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
    for res in self.protein.residues:
        for conf in res.conformers:
            for atom in conf.atoms:
                atom.connect12 = []
                for atom2 in conf.atoms:
                    if atom == atom2:
                        continue
                    if atom.distance(atom2) < 1.2:
                        atom.connect12.append(atom2)