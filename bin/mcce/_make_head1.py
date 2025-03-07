"""
Make head1.lst file for MCCE
"""

def make_head1(self, head1_file):  # self is a MCCE object
    """
    Make head1.lst file for MCCE
    head1.lst is a summary of residues and instruction for step 2.
    """
    with open(head1_file, "w") as f:
        f.write("#Rotamer Making Site Specific Instruction:\n")
        for res in self.protein.residues:
            f.write(f"{res.resname} {res.chain}{res.sequence:4d} {res.insertion} R f 00 S f 0.0 H f\n")