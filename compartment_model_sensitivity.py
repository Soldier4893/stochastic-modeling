# from reaction_network import PreNetwork, ReactionNetwork
from matplotlib import pyplot as plt
import numpy as np


class SimpleFastCompartmentManager:
    def __init__(self, ckappas, kB, kD):
        self.comparts = np.zeros(100000, dtype = np.int32)
        self.kB = np.float64(kB)
        self.kD = np.float64(kD)
        self.kI = np.float64(ckappas[0])
        self.kE = np.float64(ckappas[1])
        self.kF = np.float64(ckappas[2])
        self.kC = np.float64(ckappas[3])

    def simComparts(self, limit):        
        self.comparts[0] = 1
        self.numS = np.int32(1)
        self.population = np.zeros((100000000, 3), dtype = np.float64)
        self.numC = np.int32(1)
        self.numR = np.int32(-1)

        timer = np.float64(0)
        time_limit = np.float64(limit)

        while True:
            self.numR += 1
            self.population[self.numR,0] = timer
            self.population[self.numR,1] = self.numS
            self.population[self.numR,2] = self.numC

            ckinetic_rates = np.array(
                (self.kI, self.kE*self.numC, self.kF*self.numC,
                self.kB*self.numC, self.kD*self.numS)
            )
            rate = np.sum(ckinetic_rates)
            which_action = np.random.choice(np.arange(5), p=ckinetic_rates/rate)
            
            delay = np.random.exponential(1/rate)
            timer += delay
            if timer >= time_limit:
                break

            if which_action == 0: # add compartment
                self.add_compart(np.int32(np.random.randint(0, 11)))
            elif which_action == 1: # remove compartment
                self.remove_compart(self.random_sample_index())
            elif which_action == 2: # split compartment
                which_compart = np.random.choice(np.arange(self.numC)) # currently choice of compartment is uniform
                #which_compart = self.random_sample_index() # maybe different name
                amount = 0
                if self.comparts[which_compart] != 0:
                    amount = np.int32(np.random.randint(0, self.comparts[which_compart]+1))
                self.comparts[which_compart] -= amount
                self.numS -= amount
                self.add_compart(amount)
            elif which_action == 3: # increment one compartment's chemicals
                kinetic_rates = self.comparts[:self.numC]
                which_compart = np.random.choice(np.arange(self.numC))
                self.comparts[which_compart] += 1
                self.numS += 1
            else: # decrement one compartment's chemicals
                kinetic_rates = self.comparts[:self.numC]
                which_compart = np.random.choice(np.arange(self.numC), p = kinetic_rates/np.sum(kinetic_rates))
                self.comparts[which_compart] -= 1
                self.numS -= 1

    def getStats(self):
        return self.numC, self.numR, self.numS, self.population[:self.numR]
    
    def add_compart(self, amount):
        self.comparts[self.numC] = amount
        self.numC += 1
        self.numS += amount
    
    def remove_compart(self, index):   
        self.numS -= self.comparts[index]     
        if index != self.numC-1:
            self.comparts[index] = self.comparts[self.numC-1]            
        self.numC -= 1

    def random_sample_index(self):
        return np.int32(np.random.randint(0, self.numC))

    def graph(self):
        t = self.population[:self.numR, 0]
        s = self.population[:self.numR, 1]
        c = self.population[:self.numR, 2]
        print(self.numR)
        plt.plot(t, s)
        plt.plot(t, c)
        plt.legend(['numS','numC'])
        # plt.savefig('/foo.png', bbox_inches='tight')
        plt.show()

class CoupledCompartments:
    def __init__(self, ckappas, kB, kD):
        self.comparts = np.zeros((2, 1024**2), dtype = np.int32)
        self.kB = np.float64(kB)
        self.kD = np.float64(kD)
        self.kI = np.float64(ckappas[0])
        self.kE = np.float64(ckappas[1])
        self.kF = np.float64(ckappas[2])
        self.kC = np.float64(ckappas[3])

    def simComparts(self, time_limit, h):
        self.comparts[0][0] = 1
        self.comparts[1][0] = 1
        self.numCs = [1, 1]
        self.numSs = [1, 1]
        timer = 0

        Pks = np.random.exponential(scale=1.0, size=(6, 3))
        Tks = np.zeros((6, 3))
        Aks = np.zeros((6, 3))
        while True:
            kinetic_rates1 = np.array(
                (self.kI, self.kE*self.numCs[0], (self.kF-h)*self.numCs[0],
                  self.kC*self.numCs[0]*(self.numCs[0]-1)/2,
                    self.kB*self.numCs[0], self.kD*self.numSs[0])
            )
            kinetic_rates2 = np.array(
                (self.kI, self.kE*self.numCs[1], (self.kF+h)*self.numCs[1],
                  self.kC*self.numCs[1]*(self.numCs[1]-1)/2,
                    self.kB*self.numCs[1], self.kD*self.numSs[1])
            )

            for i in range(6):
                Aks[i][0] = min(kinetic_rates1[i], kinetic_rates2[i])
                Aks[i][1] = kinetic_rates1[i] - Aks[i][0]
                Aks[i][2] = kinetic_rates2[i] - Aks[i][0]
            
            
            with np.errstate(divide='ignore', invalid='ignore'):
                delta_tks = np.where(Aks != 0, np.divide(Pks - Tks,Aks), np.inf)
            
            mu = np.unravel_index(np.nanargmin(delta_tks, axis=None), Aks.shape)
            Tks += Aks * delta_tks[mu]
            Pks[mu] += np.random.exponential(1)
            timer += delta_tks[mu]

            if timer >= time_limit:
                break
            
            which_action = mu[0]

            if mu[1] == 0:
                self.update_comparts(which_action, 0)
                self.update_comparts(which_action, 1)
            elif mu[1] == 1:
                self.update_comparts(which_action, 0)
            elif mu[1] == 2:
                self.update_comparts(which_action, 1)

    def update_comparts(self, which_action, i):
        if which_action == 0: # add compartment 0-10
            self.add_compart(np.int32(np.random.randint(0, 11)), i)
        elif which_action == 1: # remove compartment
            self.remove_compart(self.random_sample_index(i), i)
        elif which_action == 2: # split compartment
            # randomly select compartment and record a random amount
            amount = 0
            which_compart = np.random.choice(np.arange(self.numCs[i]))
            if self.comparts[i][which_compart] != 0:
                amount = np.int32(np.random.randint(0, self.comparts[i][which_compart]+1))
            
            # create a new compartment with the amount taken away
            self.comparts[i][which_compart] -= amount
            self.numSs[i] -= amount
            self.add_compart(amount, i)
        elif which_action == 3: # merge compartments
            # randomly select compartment and record amount
            which_compart = np.random.choice(np.arange(self.numCs[i]))
            amount = self.comparts[i][which_compart]
            
            # add amount since remoce_compart automatically takes it away
            self.numSs[i] += amount
            self.remove_compart(which_compart, i)

            # select new random compartment and increment
            which_compart = np.random.choice(np.arange(self.numCs[i]))
            self.comparts[i][which_compart] += amount
        elif which_action == 4: # increment one compartment's chemicals
            which_compart = np.random.choice(np.arange(self.numCs[i]))
            self.comparts[i][which_compart] += 1
            self.numSs[i] += 1
        elif which_action == 5: # decrement one compartment's chemicals
            comparts = self.comparts[i][:self.numCs[i]]
            which_compart = np.random.choice(np.arange(self.numCs[i]), p = comparts/np.sum(comparts))
            self.comparts[i][which_compart] -= 1
            self.numSs[i] -= 1
    
    def getStats(self):
        return self.numC, self.numR, self.numS, self.population[:self.numR]
    
    def add_compart(self, amount, i):
        self.comparts[i][self.numCs[i]] = amount
        self.numCs[i] += 1
        self.numSs[i] += amount
    
    def remove_compart(self, index, i):
        self.numSs[i] -= self.comparts[i][index]     
        if index != self.numCs[i]-1:
            self.comparts[i][index] = self.comparts[i][self.numCs[i]-1]            
        self.numCs[i] -= 1

    def random_sample_index(self, i):
        return np.int32(np.random.randint(0, self.numCs[i]))
