"""
Genetic Algorithm for MCCE
"""

import logging
import random
from .constants import *

def ga_optimize(self):  # Here self is an instance of MCCE
    """
    Genetic Algorithm for MCCE
    """
    # gather constants and parameters for GA
    ga_phs = [float(self.prm.GA_PH1.value), float(self.prm.GA_PH2.value), float(self.prm.GA_PH3.value)]
    ga_pool = int(self.prm.GA_POOL.value)
    ga_maxgen = int(self.prm.GA_MAXGEN.value)
    ga_preserve = GA_preserve
    ga_single_cross = GA_single_cross
    ga_double_cross = GA_double_cross
    ga_mutation = GA_mutation
    logging.info(f"      GA pool size:{ga_pool}")
    logging.info(f"      GA maximum generations:{ga_maxgen}")
    logging.info(f"      GA solution pHs:{ga_phs}")
    logging.info(f"      GA preserve population:{ga_preserve}")
    logging.info(f"      GA single_cross rate:{ga_single_cross}")
    logging.info(f"      GA double_cross rate:{ga_double_cross}")
    logging.info(f"      GA mutation rate:{ga_mutation}")

    


