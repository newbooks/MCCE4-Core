"""
Assign charge and radius to atoms in protein object
The input pdb file doesn't have charge and radius information. This function assigns charge and radius to atoms in the protein object.
The exception is step3.py, which should be given a chance to use boundary rasii and charges from step2_out.pdb.
"""

import logging
from ._sas import radius, UNASSIGEDN_RAD

def assign_qr(self):  # Here, self is a MCCE object
    for res in self.protein.residues:
        for conf in res.conformers:
            for atom in conf.atoms:
                key = ("RADIUS", conf.conftype, atom.atomname)
                if key in self.tpl:
                    r_param = self.tpl[key]
                    atom.r_vdw = r_param.r_vdw
                    atom.e_vdw = r_param.e_vdw
                    atom.r_boundary = r_param.r_bound
                else:
                    atom.r_vdw = radius.get(atom.element, UNASSIGEDN_RAD)
                    atom.e_vdw = 0.0
                    atom.r_boundary = 2.0
                    logging.debug(f"   Radius for {key} not found, assuming r_vdw={atom.r_vdw:.2f}, e_vdw=0.00, and r_bound=2.00")
                   
                key = ("CHARGE", conf.conftype, atom.atomname)
                if key in self.tpl:
                    atom.charge = self.tpl[key]
                else:
                    atom.crg = 0.0
                    logging.debug(f"   Charge for {key} not found, assuming charge=0.00")