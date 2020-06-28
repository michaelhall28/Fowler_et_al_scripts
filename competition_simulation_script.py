import numpy as np
from clone_competition_simulation.parameters import Parameters
import pandas as pd
from clone_competition_simulation.fitness_classes import Gene, ExponentialDist, MutationGenerator
from clone_competition_simulation.sim_sampling import get_vafs_for_all_biopsies

TIMES = [2000, 3000, 4000]
DIVISION_RATE = 0.5
COVERAGE = 690  # Average coverage.
DETECTION_LIMIT = 5

BIOPSY_SHAPE = 250, 120  # 30000 basal cells per biopsy.

biopsy_grid = 8, 5  # 40 biopsies per "patient"

NUM_SIM_SAMPLES = biopsy_grid[0]*biopsy_grid[1]
GRID_SHAPE=BIOPSY_SHAPE[0]*biopsy_grid[0], BIOPSY_SHAPE[1]*biopsy_grid[1]
TOTAL_CELLS = GRID_SHAPE[0]*GRID_SHAPE[1]
BIOPSIES = []
for i in range(biopsy_grid[0]):
    for j in range(biopsy_grid[1]):
        BIOPSIES.append({'biopsy_origin': (i*BIOPSY_SHAPE[0], j*BIOPSY_SHAPE[1]), 'biopsy_shape': BIOPSY_SHAPE})

def run_sim(mutation_generator, mutation_rate):
    p = Parameters(algorithm='WF2D', grid_shape=GRID_SHAPE,
                   mutation_generator=mutation_generator,
                   mutation_rates=mutation_rate, times=TIMES,
                   print_warnings=False, division_rate=DIVISION_RATE, samples=2,
                   cell_in_own_neighbourhood=True)
    s = p.get_simulator()
    s.run_sim()
    res = {t: get_vafs_for_all_biopsies(s, BIOPSIES, COVERAGE, DETECTION_LIMIT, merge_clones=False, sample_num=i) for i, t in enumerate(TIMES)}
    return res


if __name__ == "__main__":
    # So that simulations can be run in parallel
    # this runs a set of simulations with the same mutation rate and driver proportions
    # The argument to the script is simply a number between 1 and 30
    # which determines which parameter combination to run
    import sys
    sim_num = sys.argv[1]
    sim_num = int(sim_num)
    mutation_rates = [0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01]
    driver_proportions = [0.00001, 0.0001, 0.001, 0.01, 0.1]
    
    driver_fitnesses = [1.01, 1.05, 1.1, 1.2, 1.5, 2, 5]
    
    count = 0
    for m in mutation_rates:
        for p in driver_proportions:
            count += 1
            if count == sim_num:
                for f in driver_fitnesses:
                    np.random.seed(0)
                    mutation_generator = MutationGenerator(genes=[Gene('driver', ExponentialDist(mean=f, offset=1),  # Fitness drawn from offset + Exp(mean-offset)
                                                           synonymous_proportion=1-p)],
                                                           combine_mutations='multiply', 
                                                           multi_gene_array=False)

                    clones = run_sim(mutation_generator, m)
                    for t, c in clones.items():
                        c.to_csv("clones_exp_fitness_mult_{}_{}_{}_{}.csv".format(m, p, f, t))


                          
