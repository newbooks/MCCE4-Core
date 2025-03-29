"""
GA optimization for MCCE
This is the main module imported by main.py that developers can implement with their own strategy. 
It would use the following provided modules:
_ga_lib.py: utility functions for GA (not encourged to modify)
_ga_forcefield.py: forcefield functions for GA (not encourged to modify)
"""

import logging
from .constants import *
from ._ga_lib import *

def ga_optimize(self, writepdb=False):  # Here self is an instance of MCCE
    """
    Genetic Algorithm for MCCE
    """
    # gather constants and parameters for GA
    ga_pool_size = int(self.prm.GA_POOL.value)
    ga_max_generations = int(self.prm.GA_MAXGEN.value)
    # the following two variables are numbers of crossovers and mutations in one period (generation)
    ga_crossover_number = int(GA_crossover_rate * ga_pool_size)
    ga_mutation = int(GA_mutation_rate * ga_pool_size)

    # Save the current handlers
    logger = logging.getLogger()
    original_handlers = logger.handlers.copy()
    # Create a file handler
    file_handler = logging.FileHandler(GA_PROGRESS, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s", datefmt='%Y-%m-%d %H:%M:%S'))
    logger.handlers = [file_handler]
    # From now on, the log goes to the file GA_PROGRESS


    # Start the GA optimization
    logging.info("Genetic Algorithm optimization for MCCE:")
    logging.info(f"   GA pool size:{ga_pool_size}")
    logging.info(f"   GA maximum generations:{ga_max_generations}")
    logging.info(f"   GA crossover rate:{GA_crossover_rate}")
    logging.info(f"   GA mutation rate:{GA_mutation_rate}")
    logging.info(f"   Prepare GA pool ...")

    # Initialize the pool, this includes computing fitness score for each individual and sorting the pool
    pool = Pool(mcce=self, size=ga_pool_size)
    # pool.population[0].print_connect12()
    logging.info(f"      Done preparing pool.")

    # testing the pool clone function and time
    logging.info(f"      Start to clone testing ...")
    test_clone(pool)            
    logging.info(f"      Done cloning and testing.")

    # Start the evolution
        # Crossover, loop over ga_crossover_number times
            # select two parents
            # crossover -  obtain two new individuals
            # calculate fitness scores for the new individuals
            # replace the worst individuals with the new individuals
        # Mutation, loop over ga_mutation times
            # select one individual
            # mutate - obtain a new individual
            # calculate fitness score for the new individual
            # replace the worst individual with the new individual
        # Compute PFA and test termination condition
    # End of evolution loop, we get an updated pool with the best individuals



    # logging.info(f"      Start to evolve ...")
    # pool.evolve()
    # logging.info(f"      Done evolution.")
    
    
    
    if writepdb:
        logging.info(f"      Write individuals to {GA_OUTPUT_FOLDER}")
        pool.writepdb()
    logging.info("Genetic Algorithm optimization completed.")

    # pool.population[2].test_lineage()
    
    
    # Restore the original logger handlers
    logger.handlers = original_handlers