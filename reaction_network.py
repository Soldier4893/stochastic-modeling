# modified from GPT-3.5 generated code
import numpy as np
import scipy
import prenetwork
from matplotlib import pyplot as plt

class reactionNetwork:
    def __init__(self, X, zeta_table, rate_table, rates, time_steps):
        self.X = X
        self.k = len(rates)
        self.n = len(rate_table[0])
        self.zeta_tables = zeta_table
        self.rate_table = rate_table
        self.rates = rates
        self.time_steps = time_steps

    def simulate(self):
        for _ in range(self.time_steps):
            rate = self.calculate_total_rate()
            delay = np.random.exponential(1.0/rate)

    def calculate_total_rate(self):
        kinetics = np.zeros(len(rates))
        for kinetic_rate in rate_table:
            # iterate over n as up to n chemicals may be involved in the reaction
            for i in range(self.n):
                kinetics[i] = self.kinetic_function(i, kinetic_rate[i])

    def kinetic_function(self, i, kinetic_rate):
        if kinetic_rate == 0:
            return 0
        return scipy.special.comb(self.X[i], kinetic_rate)


# Example usage:
reactions = [
    ("2A+B","C+D", 0.9),
    ("C","E", 0.5),
    ("D+E","F", 10)
]

species = ['A','B','C','D','E','F']
X = np.array([10, 20, 0, 0, 0, 0])
assert len(species) == len(X), "species set not aligned with initial population"
time_steps = 5
pN = prenetwork.preNetwork(species, reactions)
zeta_table, rate_table, rates = pN.get_tables()
print(zeta_table, '\n\n',rate_table, '\n\n', rates)



# network = reactionNetwork(X, zeta_table, rate_table, rates, time_step)
# network.simulate(initial_concentrations)
