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
        self.population_table = np.zeros((time_steps, self.n))
        self.times = np.zeros(time_steps+1)

    def simulate(self):
        for i in range(self.time_steps):
            # record population at this time
            self.population_table[i] = self.X
            
            # calculate kinetic rates
            kinetic_rates = self.calculate_kinetics() * self.rates
            rate = np.sum(kinetic_rates)

            # calculate delay and time variable
            delay = np.random.exponential(1.0/rate)
            self.times[i+1] = self.times[i] + delay

            # pick reaction and add to X
            which_reaction = np.random.choice(np.arange(self.k), p=kinetic_rates/rate)
            self.X += self.zeta_table[which_reaction, :]

    def calculate_kinetics(self):
        kinetics = np.ones(self.k)
        for i in range(self.n):
            # iterate over n as up to n species may be involved in the reaction
            for species_index in range(self.n):
                kinetics[i] *= self.kinetic_function(species_index, rate_table[i, species_index])
        return kinetics

    def kinetic_function(self, species_index, kinetic_rate):
        if kinetic_rate == 0:
            return 0
        return scipy.special.comb(self.X[species_index], kinetic_rate)
    
    def graph(self):
        for i in range(self.n):
            plt.plot(self.times, self.population_table[:, i])

# Example usage:
reactions = [
    ("2A+B","C+D", 0.9),
    ("C","E", 0.5),
    ("D+E","F", 10)
]

species = ['A','B','C','D','E','F']
X = np.array([10, 20, 0, 0, 0, 0], dtype = np.uint16)
assert len(species) == len(X), "species set not aligned with initial population"
time_steps = 5
pN = prenetwork.preNetwork(species, reactions)
zeta_table, rate_table, rates = pN.get_tables()
print(zeta_table[0], '\n\n',rate_table, '\n\n', rates)



# network = reactionNetwork(X, zeta_table, rate_table, rates, time_step)
# network.simulate(initial_concentrations)
