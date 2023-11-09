from compartment_manager import CompartmentManager
from reaction_network import PreNetwork
import numpy as np

# reactions = [
#     ("G","G+M", 200),
#     ("M","M+P", 10),
#     ("M","0", 25),
#     ("P","0", 1),
#     ("2P","D", 0.01),
#     ("D","0", 1)
# ]
reactions = [("0", "S", 1), ("S", "0", 1)]
species = ['S', 'B']
X = np.array([10, 9], dtype = np.int32)
ctime_steps = 1500

def sample_pmf_X():
    # temporary distribution
    return np.array((np.int32(np.random.randint(0, 11)), 0))
    # return np.array([1, np.random.randint(0, 11), np.random.randint(0, 11), np.random.randint(0, 11)], dtype = np.int32)

if __name__ == "__main__":
    PN = PreNetwork(species, reactions)
    ckappas = np.array([1, 1, 2, 0]) #in order is I, E, F, C
    assert len(species) == len(X), "species set not aligned with initial population"
    
    cm = CompartmentManager(len(species), PN, sample_pmf_X, ctime_steps, ckappas)
    cm.simCompartments()
    cm.graph()


