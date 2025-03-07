"""
Load mccepdb file from step1_out.pdb
"""

from .constants import *
from .pdbio import *


def load_mccepdb(self):  # self is an instance of MCCE
    atoms = []
    with open(STEP1_OUT, "r") as f:
        # read line by line from f
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = Atom()
                atom.load_mcceline(line)
                atoms.append(atom)
    
    protein = Protein()
    for atom in atoms:
        residue = None
        for res in protein.residues:
            if res.resid == (atom.resname, atom.chain, atom.sequence, atom.insertion):
                residue = res
                break
        if residue is None:
            residue = Residue()
            residue.resname = atom.resname
            residue.chain = atom.chain
            residue.sequence = atom.sequence
            residue.insertion = atom.insertion
            residue.resid = (atom.resname, atom.chain, atom.sequence, atom.insertion)
            protein.residues.append(residue)

        # create a new conformer if the conformer does not exist
        conformer = None
        for conf in residue.conformers:
            if conf.confid == f"{atom.resname}{atom.chain}{atom.sequence:04d}{atom.insertion}{atom.confnum:03d}":
                conformer = conf
                break
        if conformer is None:
            conformer = Conformer()
            conformer.confnum = atom.confnum
            conformer.confid = f"{atom.resname}{atom.chain}{atom.sequence:04d}{atom.insertion}{atom.confnum:03d}" 
            conformer.conftype = atom.resname + atom.history[:2]
            conformer.history = atom.history
            conformer.altloc = atom.altloc
            conformer.resname = atom.resname
            conformer.chain = atom.chain
            conformer.sequence = atom.sequence
            conformer.insertion = atom.insertion
            conformer.parent_residue = residue
            residue.conformers.append(conformer)

        # add this atom to conformer
        atom.parent_conf = conformer
        conformer.atoms.append(atom)

    # lastly, scan all residues to see if the first conformer is the backbone
    for res in protein.residues:
        if res.conformers[0].conftype[-2:] != "BK":
            bk_conf = Conformer()
            bk_conf.confnum = 0
            bk_conf.confid = f"{res.resname}{res.chain}{res.sequence:04d}{res.insertion}000"
            bk_conf.conftype = res.resname + "BK"
            bk_conf.history = "BKO000_000"
            bk_conf.parent_residue = res
            bk_conf.atoms = []
            res.conformers.insert(0, bk_conf)

    self.protein = protein