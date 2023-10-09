from reaction_network import PreNetwork, ReactionNetwork
from matplotlib import pyplot as plt
import numpy as np

np.random.seed(18)

class CompartmentManager:
    def __init__(self, n, pN, mu_X, ctime_steps, ckappas):
        self.d = {}
        self.compartments = []
        self.pN = pN
        self.mu_X = mu_X
        self.ckappas = ckappas
        self.ctimes = np.zeros(ctime_steps)
        self.ctime_steps = ctime_steps
        self.cpopulation = np.zeros(ctime_steps, dtype = np.int64)
        self.n = n
    
    def simCompartments(self):
        zeta_table, kappa_table, kappas = self.pN.get_tables()
        self.populations = np.zeros((self.ctime_steps, len(kappa_table[0])))

        # initial amount of compartments
        for _ in range(3):
            self.add(ReactionNetwork(self.mu_X(), zeta_table, kappa_table, kappas))
        self.cpopulation[0] = 3

        for i in range(1, self.ctime_steps):
            nc = len(self.compartments)
            kinetic_rates = np.array([1, nc, nc, nc*(nc-1)/2]) * self.ckappas
            if not kinetic_rates.any():
                print("no more reactions")
                break

            rate = np.sum(kinetic_rates)
            which_action = np.random.choice(np.arange(4), p=kinetic_rates/rate)

            delay = np.random.exponential(1/rate)
            self.ctimes[i] = self.ctimes[i-1] + delay
            self.cpopulation[i] = len(self.compartments)

            if which_action == 0: # add compartment
                self.add(ReactionNetwork(self.mu_X(), zeta_table, kappa_table, kappas))
            elif which_action == 1: # remove compartment
                self.pop()
            elif which_action == 2: # split compartment
                cpmt_to_split = self.sample()
                self.add(ReactionNetwork(cpmt_to_split.split_X(), zeta_table, kappa_table, kappas))
            else: #which_action == 3  # merge compartments
                cpmt_giver = self.pop()
                cpmt_receiver = self.sample()
                cpmt_receiver.merge_X(cpmt_giver.X)

            for compartment in self.compartments:
                compartment.simGillespie(delay)
                self.populations[i] += compartment.X

    def graph(self):
        for i in range(self.n):
            plt.plot(self.ctimes, self.populations[:, i])
        plt.plot(self.ctimes, self.cpopulation)
        plt.legend(['G','M','P','D','C']) # change later
        plt.show()

    def add(self, element):
        # enforce unique items
        # if self.d[element]:
        #     return
        self.compartments.append(element)
        self.d[element] = len(self.compartments)-1
    
    def remove_and_return(self, element):   
        i = self.d[element]
        del self.d[element]
        if i == len(self.compartments)-1:
            self.compartments.pop()
        else:
            self.compartments[i] = self.compartments.pop()
            self.d[self.compartments[i]] = i
        return element

    def sample(self):
        return np.random.choice(self.compartments) #p = self.generate_distribution

    # normalize before returning
    def generate_distribution(self):
        return np.ones(len(self.compartments), dtype = np.float64)/len(self.compartments) #perhaps change later

    def pop(self):
        temp = self.sample()
        return self.remove_and_return(temp)
    
    