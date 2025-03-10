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
    self.reset_connect12()
    for i_res in range(len(self.protein.residues)): # Loop over residues and we need index number to look up neighboring residues
        res = self.protein.residues[i_res]
        for conf in res.conformers:
            for atom in conf.atoms:
                key = ("CONNECT", atom.atomname, conf.conftype)
                connected_atomnames = self.tpl.get(key, [])
                print(f"CONNECT {atom.atomname} {conf.conftype} {connected_atomnames}")