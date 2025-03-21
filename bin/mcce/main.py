"""
Define MCCE object and its methods
This is the main module of MCCE. It contains the MCCE class and its methods.
"""

from ._rot_stat import RotStat

class MCCE:
    def __init__(self, prm=None, tpl=None, protein=None):
        self.prm = prm
        self.tpl = tpl
        self.protein = protein
        self.rotate_rules = {}

    from ._strip_cofactors import strip_cofactors
    from ._sas import sas_protein
    from ._assign_qr import assign_qr
    from ._make_head1 import make_head1
    from ._load_mccepdb import load_mccepdb
    from ._place_missing_heavy_atoms import place_missing_heavy_atoms
    from ._make_connect import reset_connect12
    from ._make_connect import make_connect12
    from ._make_connect import make_connect12_fast
    from ._make_connect import print_connect12
    from ._make_connect import check_connect12
    from ._make_connect import make_connect13
    from ._make_connect import print_connect13
    from ._make_connect import check_connect13
    from ._make_connect import make_connect14
    from ._make_connect import check_connect14
    from ._rotate import prepare_rotate_rules
    from ._rotate import print_rotate_rules
    from ._rotate import apply_rotate_rules
    from ._conftype_conformers import propogate_swap
    from ._conftype_conformers import propogate_conftypes
    from ._ga import ga_optimize