"""
Define MCCE object and its methods
This is the main module of MCCE. It contains the MCCE class and its methods.
"""

class MCCE:
    def __init__(self, prm=None, tpl=None, protein=None):
        self.prm = prm
        self.tpl = tpl
        self.protein = protein

    from ._strip_cofactors import strip_cofactors