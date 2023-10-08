from compartment_manager import CompartmentManager
from reaction_network import PreNetwork
import numpy as np

reactions = [
    ("G","G+M", 200),
    ("M","M+P", 10),
    ("M","0", 25),
    ("P","0", 1),
    ("2P","D", 0.01),
    ("D","0", 1)
]
species = ['G','M','P', 'D']
X = np.array([1, 0, 0, 0], dtype = np.int64)
ctime_steps = 200

def sample_pmf_X():
    # temporary distribution
    return np.array([1, np.random.randint(0, 11), np.random.randint(0, 11), np.random.randint(0, 11)], dtype = np.int32)

if __name__ == "__main__":
    PN = PreNetwork(species, reactions)
    ckappas = np.array([0.2, 0.3, 0.4, 0.5]) #in order is I, E, F, C
    assert len(species) == len(X), "species set not aligned with initial population"
    
    cm = CompartmentManager(PN, sample_pmf_X, ctime_steps, ckappas)
    cm.simCompartments()
    cm.graph()


