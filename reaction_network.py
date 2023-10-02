# modified from GPT-3.5 generated code
import numpy as np
import scipy
import re
from matplotlib import pyplot as plt

# given human readable reaction network output as np arrays zetas, kappa_table
# species is list, reaction is list of tuples as (reactants, products, rate) as string, string, float.
# kappa_table is n x k, kappa_table[i][j] tells how much of species j is needed for reaction i
# zeta_table is n x k, when reaction i occurs we should set X += zeta_table[i]
class preNetwork:
    def __init__(self, species, reactions):
        n, k = len(species), len(reactions)
        self.zeta_table = np.zeros((k, n), dtype = np.int32)
        self.kappa_table = np.zeros((k, n), dtype = np.int32)
        self.rates = np.zeros(k)
        self.species_dict = {item: index for index, item in enumerate(species)}

        for reaction_num in range(k):
            reactants, products, rate = reactions[reaction_num]
            self.rates[reaction_num] = rate
            
            # find reactants, add to zeta
            self.addSideToZeta(reactants, reaction_num, sign = -1)
            self.addSideToZeta(products, reaction_num, sign = 1)

    def addSideToZeta(self, side, reaction_num, sign = -1):
        zeta_indices = self.parse(side)
        for species_index, amount in zeta_indices:
            self.zeta_table[reaction_num, species_index] += amount * sign
            self.kappa_table[reaction_num, species_index] += amount * (sign == -1)

    def parse(self, string):
        zeta_indices = []
        complexes_in_side = string.replace(" ","").split("+")
        for complex in complexes_in_side:
            zeta_indices.append(self.complex_to_index(complex))
        return zeta_indices

    def complex_to_index(self, complex):
        if complex == '0':
            return [], 0
        
        numbers_as_strings = re.findall(r'\d+', complex)
        amount = 0
        if numbers_as_strings:
            amount = int(numbers_as_strings[0])
        else:
            amount = 1
        species_index = self.species_dict[re.sub(r'[^a-zA-Z]', '', complex)]
        return species_index, amount
    
    def get_tables(self):
        return self.zeta_table, self.kappa_table, self.rates

class reactionNetwork:
    def __init__(self, X, zeta_table, kappa_table, rates, time_steps):
        self.X = X
        self.k = len(rates)
        self.n = len(kappa_table[0])
        self.zeta_table = zeta_table
        self.kappa_table = kappa_table
        self.rates = rates
        self.time_steps = time_steps
        self.population_table = np.zeros((time_steps, self.n))
        self.times = np.zeros(time_steps)

    def simGillespie(self):
        # initial population
        self.population_table[0] = self.X

        for i in range(1, self.time_steps):
            # record population at this time
            self.population_table[i] = self.X
            
            # calculate kinetic rates
            kinetic_rates = self.calculate_kinetics() * self.rates
            if not kinetic_rates.any():
                print("no more reactions")
                break

            rate = np.sum(kinetic_rates)

            # calculate delay and time variable
            delay = np.random.exponential(1.0/rate)
            self.times[i] = self.times[i-1] + delay

            # pick reaction and add to X
            which_reaction = np.random.choice(np.arange(self.k), p=kinetic_rates/rate)
            self.X += self.zeta_table[which_reaction, :]
    
    def simNextReaction(self):
        # initial population
        self.population_table[0] = self.X
        T_ks = np.zeros(self.k)
        for i in range(1, self.time_steps):
            # record population at this time
            self.population_table[i] = self.X
            
            # calculate kinetic rates
            kinetic_rates = self.calculate_kinetics() * self.rates
            if not kinetic_rates.any():
                print("no more reactions")
                break

            # find delta := min delta_k
            P_ks = np.log(1/np.random.rand(self.k))         
            delta = np.nanmin((P_ks - T_ks)/kinetic_rates)
            
            # increment time, X
            self.times[i] = self.times[i-1] + delta
            self.X += self.zeta_table[9] # not 9, should be mu. Not sure what mu is.
            T_ks += kinetic_rates * delta
            P_ks[9] += np.log(1/np.random.rand())
            
    def calculate_kinetics(self):
        kinetics = np.ones(self.k)
        for i in range(self.k):
            # up to n species may be involved in the reaction
            for j in range(self.n):
                # if the species is relevant (a.k.a. self.kappa_table[i, j] != 0) then the rate is multiplied accordingly
                # otherwise rate is multiplied by one
                kinetics[i] *= scipy.special.comb(self.X[j], self.kappa_table[i, j]) * (self.kappa_table[i, j] != 0) + (self.kappa_table[i, j] == 0)
        return kinetics
    
    def graph(self):
        for i in range(self.n):
            plt.plot(self.times, self.population_table[:, i])
        plt.show()

class nextReactionSim:
    def __init__(self, X, zeta_table, kappa_table, rates, time_steps):
        pass


# # Example usage:
# reactions = [
#     ("G","G+M", 200),
#     ("M","M+P", 10),
#     ("M","0", 25),
#     ("P","0", 1),
#     ("2P","D", 0.01),
#     ("D","0", 1)
# ]

# species = ['G','M','P', 'D']
# X = np.array([1, 0, 0, 0], dtype = np.int32)
# assert len(species) == len(X), "species set not aligned with initial population"
# time_steps = 5000
# pN = preNetwork(species, reactions)
# zeta_table, kappa_table, rates = pN.get_tables()

# network = reactionNetwork(X, zeta_table, kappa_table, rates, time_steps)
# network.simulate()
# network.graph()