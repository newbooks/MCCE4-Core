"""
Place missing heavy atoms in the structure
"""

def place_missing_heavy_atoms(self):    # self is a MCCE object
    n_placed = 0

    # make 12 connectivity
    self.reset_connect12()
    self.make_connect12()



    return n_placed