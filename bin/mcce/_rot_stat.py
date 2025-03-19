"""
Step 2 Rotamer making statistics
"""

import logging
from .constants import ROT_STAT

class RotStat:
    class StatItem:
        def __init__(self):
            self.start = 0      # conformer count ath the start
            self.swap = 0       # conformer count after swap
            self.type = 0       # conformer count after propogating conformer types
            self.ga = 0         # conformer count after GA selection 80% of population
            self.prune = 0      # conformer count after pruning by intra vdw and backbone vdw
            self.cluster = 0    # conformer count after pruning and clustering
        
    def __init__(self, protein):
        self.res_rot_stat = [RotStat.StatItem() for res in protein.residues]  # initialize a list of StatItem for each residue

    def count_stat(self, protein, step=None):
        """
        Count for step named in step
        """
        attributes = [x for x in self.res_rot_stat[0].__dict__.keys()]
        if step in attributes:
            for ires in range(len(protein.residues)):
                n_conf = len(protein.residues[ires].conformers) - 1  # -1 because the first conformer is the backbone
                if n_conf < 0:
                    n_conf = 0
                    logging.warning(f"Residue {protein.residues[ires].resid} has no conformer")
                setattr(self.res_rot_stat[ires], step, n_conf)
        else:
            logging.error(f"Step {step} is not in the list of attributes: {attributes}")
    
    def write_stat(self, protein):
        header = "  Residue   Start    Swap    Type      GA   Prune Cluster\n"
        total_conf = self.StatItem()
        attributes = [x for x in self.res_rot_stat[0].__dict__.keys()]
        """
        Write rotamer statistics to a file
        """
        with open(ROT_STAT, "w") as f:
            f.write(header)
            for ires in range(len(protein.residues)):
                res = protein.residues[ires]
                stat = self.res_rot_stat[ires]
                f.write(f"{res.resname}{res.chain}{res.sequence:04d}{res.insertion:1s}")
                for attr in attributes:
                    f.write(f"{stat.__dict__[attr]:8d}")
                    total_conf.__dict__[attr] += stat.__dict__[attr]
                f.write("\n")
            f.write("-" * len(header) + "\n")
            f.write(f"{'Total':>9}")
            for attr in attributes:
                f.write(f"{total_conf.__dict__[attr]:8d}")
            f.write("\n")