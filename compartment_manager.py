from reaction_network import PreNetwork, ReactionNetwork
from matplotlib import pyplot as plt
import numpy as np

np.random.seed(19)

class CompartmentManager:
    def __init__(self, n, pN, mu_X, ctime_steps, ckappas):
        self.d = {}
        self.compartments = []
        self.pN = pN
        self.mu_X = mu_X
        self.ckappas = ckappas
        self.ctimes = np.zeros(ctime_steps)
        self.ctime_steps = ctime_steps
        self.cpopulation = np.zeros(ctime_steps, dtype = np.int32)
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
    
class SimpleFastCompartmentManager:
    def __init__(self, ckappas, kB, kD):
        self.comparts = np.zeros(20000, dtype = np.int32) #max 5000 compartments, change as needed
        self.kB = np.float64(kB)
        self.kD = np.float64(kD)
        self.kI = np.float64(ckappas[0])
        self.kE = np.float64(ckappas[1])
        self.kF = np.float64(ckappas[2])
        self.kC = np.float64(ckappas[3])

        self.comparts[0] = 5
        self.comparts[1] = 3
        self.comparts[2] = 2
        self.numcomparts = np.int32(3)
        self.numS = np.int32(10)
        self.numzeros = np.int32(0)
    
    def simComparts(self, limit):
        timer = np.float64(0)
        time_limit = np.float64(limit)
        # initial amount of compartments
        self.numcomparts = np.int32(3)

        while timer < time_limit:
            ckinetic_rates = np.array(
                (self.kI, self.numcomparts*self.kE, self.kF*np.sum(self.comparts[0:self.numcomparts]),
                self.numcomparts*self.kB, self.kD*np.sum(self.comparts[:self.numcomparts])) #sum(comparts[:numcomparts+1]) is total number of S chemicals
            )
            rate = np.sum(ckinetic_rates)
            which_action = np.random.choice(np.arange(5), p=ckinetic_rates/rate)
            delay = np.random.exponential(1/rate)
            timer += delay

            if which_action == 0: # add compartment
                self.add_compart(np.int32(np.random.randint(0, 11)))
            elif which_action == 1: # remove compartment
                self.remove_compart(self.random_sample_index())
            elif which_action == 2: # split compartment
                which_compart = np.random.choice(np.arange(self.numcomparts), p = self.comparts[:self.numcomparts]/np.sum(self.comparts[0:self.numcomparts]))
                #which_compart = self.random_sample_index() # maybe different name
                amount = 0
                if self.comparts[which_compart] != 0:
                    amount = np.int32(np.random.randint(0, self.comparts[which_compart]))
                self.comparts[which_compart] -= amount
                self.add_compart(amount)
            elif which_action == 3: # increment one compartment's chemicals
                kinetic_rates = self.comparts[:self.numcomparts]
                which_compart = np.random.choice(np.arange(self.numcomparts))
                self.comparts[which_compart] += 1
            else: # decrement one compartment's chemicals
                kinetic_rates = self.comparts[:self.numcomparts]
                which_compart = np.random.choice(np.arange(self.numcomparts), p = kinetic_rates/np.sum(kinetic_rates))
                self.comparts[which_compart] -= 1
        return self.comparts[:self.numcomparts]
    
    def add_compart(self, amount):
        self.comparts[self.numcomparts] = amount
        self.numcomparts += 1
    
    def remove_compart(self, index):        
        if index != self.numcomparts-1:
            self.comparts[index] = self.comparts[self.numcomparts-1]            
        self.numcomparts -= 1

    def random_sample_index(self):
        return np.int32(np.random.randint(0, self.numcomparts))

    def graph(self):
        for i in range(self.n):
            plt.plot(self.ctimes, self.populations[:, i])
        plt.plot(self.ctimes, self.cpopulation)
        plt.legend(['G','M','P','D','C']) # change later
        plt.show()

class SimpleFastCompartmentManager1:
    def __init__(self, ckappas, kB, kD):
        self.comparts = np.zeros(20000, dtype = np.int32) #max 5000 compartments, change as needed
        self.kB = np.float64(kB)
        self.kD = np.float64(kD)
        self.kI = np.float64(ckappas[0])
        self.kE = np.float64(ckappas[1])
        self.kF = np.float64(ckappas[2])
        self.kC = np.float64(ckappas[3])

        self.comparts[0] = 5
        self.comparts[1] = 3
        self.comparts[2] = 2
        self.numcomparts = np.int32(3)
        self.numS = np.int32(10)
        self.numzeros = np.int32(0)
    
    def simComparts(self, limit):
        timer = np.float64(0)
        time_limit = np.float64(limit)
        # initial amount of compartments
        self.numcomparts = np.int32(3)
        o=0
        while timer < time_limit:
            ckinetic_rates = np.array(
                (self.kI, self.numcomparts*self.kE, self.kF*self.numS,
                self.numcomparts*self.kB, self.kD*self.numS) #sum(comparts[:numcomparts]) is total number of S chemicals
            )
            rate = np.sum(ckinetic_rates)
            which_action = np.random.choice(np.arange(5), p=ckinetic_rates/rate)
            
            delay = np.random.exponential(1/rate)
            timer += delay

            if which_action == 0: # add compartment
                self.add_compart(np.int32(np.random.randint(0, 11)))
            elif which_action == 1: # remove compartment
                self.remove_compart(self.random_sample_index())
            elif which_action == 2: # split compartment
                which_compart = np.random.choice(np.arange(self.numcomparts), p = self.comparts[:self.numcomparts]/np.sum(self.comparts[0:self.numcomparts]))
                #which_compart = self.random_sample_index() # maybe different name
                amount = 0
                if self.comparts[which_compart] != 0:
                    amount = np.int32(np.random.randint(0, self.comparts[which_compart]))
                self.comparts[which_compart] -= amount
                self.numS -= amount
                self.add_compart(amount)
            elif which_action == 3: # increment one compartment's chemicals
                kinetic_rates = self.comparts[:self.numcomparts]
                which_compart = np.random.choice(np.arange(self.numcomparts))
                self.comparts[which_compart] += 1
                self.numS += 1
            else: # decrement one compartment's chemicals
                kinetic_rates = self.comparts[:self.numcomparts]
                which_compart = np.random.choice(np.arange(self.numcomparts), p = kinetic_rates/np.sum(kinetic_rates))
                self.comparts[which_compart] -= 1
                self.numS -= 1
                o+=1
        # print(ckinetic_rates, self.numS)

        return self.comparts[:self.numcomparts],o
    
    def add_compart(self, amount):
        self.comparts[self.numcomparts] = amount
        self.numcomparts += 1
        self.numS += amount
    
    def remove_compart(self, index):   
        self.numS -= self.comparts[index]     
        if index != self.numcomparts-1:
            self.comparts[index] = self.comparts[self.numcomparts-1]            
        self.numcomparts -= 1

    def random_sample_index(self):
        return np.int32(np.random.randint(0, self.numcomparts))

    def graph(self):
        for i in range(self.n):
            plt.plot(self.ctimes, self.populations[:, i])
        plt.plot(self.ctimes, self.cpopulation)
        plt.legend(['G','M','P','D','C']) # change later
        plt.show()

class SimpleFastCompartmentManager2:
    def __init__(self, ckappas, kB, kD):
        self.comparts = np.zeros(20000, dtype = np.int32) #max 5000 compartments, change as needed
        self.kB = np.float64(kB)
        self.kD = np.float64(kD)
        self.kI = np.float64(ckappas[0])
        self.kE = np.float64(ckappas[1])
        self.kF = np.float64(ckappas[2])
        self.kC = np.float64(ckappas[3])

        self.comparts[0] = 5
        self.comparts[1] = 3
        self.comparts[2] = 2
        self.numcomparts = np.int32(3)
        self.numS = np.int32(10)
        self.numzeros = np.int32(0)
    
    def simComparts(self, limit):
        timer = np.float64(0)
        time_limit = np.float64(limit)
        # initial amount of compartments
        o=0
        while timer < time_limit:
            ckinetic_rates = np.array(
                (self.kI, (self.numcomparts+self.numzeros)*self.kE, self.kF*self.numS,
                (self.numcomparts+self.numzeros)*self.kB, self.kD*self.numS)
            )
            rate = np.sum(ckinetic_rates)
            which_action = np.random.choice(np.arange(5), p=ckinetic_rates/rate)
            delay = np.random.exponential(1/rate)
            timer += delay

            if which_action == 0: # add compartment
                self.add_compart(np.int32(np.random.randint(0, 11)))

            elif which_action == 1: # remove compartment
                which_compart = np.int32(np.random.randint(0, self.numcomparts+self.numzeros))
                self.remove_compart(which_compart)

            elif which_action == 2: # split compartment
                which_compart = np.random.choice(np.arange(self.numcomparts), p = self.comparts[:self.numcomparts]/self.numS)
                amount = 0
                if self.comparts[which_compart] != 0:
                    amount = np.int32(np.random.randint(0, self.comparts[which_compart]))
                self.comparts[which_compart] -= amount
                self.numS -= amount
                self.add_compart(amount)

            elif which_action == 3: # increment one compartment's chemicals
                kinetic_rates = self.comparts[:self.numcomparts]
                which_compart = np.random.choice(np.arange(self.numcomparts+self.numzeros))
                
                # selected empty compartment
                if which_compart >= self.numcomparts:
                    self.numzeros -= 1
                    self.add_compart(1)
                
                # selected normal compartment
                else:
                    self.comparts[which_compart] += 1
                    self.numS += 1

            else: # decrement one compartment's chemicals
                kinetic_rates = self.comparts[:self.numcomparts]
                which_compart = np.random.choice(np.arange(self.numcomparts), p = kinetic_rates/np.sum(kinetic_rates))
                self.comparts[which_compart] -= 1
                self.numS -= 1
                if self.comparts[which_compart] == 0:
                    self.remove_compart(which_compart)
            o += 1
        # print("here", self.comparts[:self.numcomparts], self.numS)
        return self.comparts[:self.numcomparts], self.numzeros, o
    
    def add_compart(self, amount):
        # add fake compartment
        if amount == 0:
            self.numzeros += 1
        # add real compartment
        else:
            self.comparts[self.numcomparts] = amount
            self.numcomparts += 1
            self.numS += amount
    
    def remove_compart(self, index):      
        # remove zero
        if index >= self.numcomparts:
            self.numzeros -= 1
            return
        
        # remove S
        self.numS -= self.comparts[index]
        
        # replace compartment
        if index < self.numcomparts-1:
            self.comparts[index] = self.comparts[self.numcomparts-1]

        # remove compartment
        self.numcomparts -= 1

    def graph(self):
        for i in range(self.n):
            plt.plot(self.ctimes, self.populations[:, i])
        plt.plot(self.ctimes, self.cpopulation)
        plt.legend(['G','M','P','D','C']) # change later
        plt.show()
