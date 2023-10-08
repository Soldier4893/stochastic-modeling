from reaction_network import PreNetwork, ReactionNetwork
from compartment_manager import CompartmentManager
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
X = np.array([1, 0, 0, 0], dtype = np.int32)
ctime_steps = 500

def sample_pmf_X():
    # temporary distribution
    return np.array([1, np.random.randint(0, 11), np.random.randint(0, 11), np.random.randint(0, 11)])

def simulate_compartments(pN, mu_X, num_compartments, ctime_steps, ckappas):
    zeta_table, rate_table, rates = pN.get_tables()

    # initial amount of compartments
    for i in range(num_compartments):
        network = ReactionNetwork(mu_X(), zeta_table, rate_table, rates)
    
    
    for _ in range(ctime_steps):
        which_action = np.random.choice(np.arange(4), p=ckappas/np.sum(ckappas))
        which_compartment = np.random.choice
        if which_action == 0:
            pass
        elif which_action == 1:
            pass
        elif which_action == 2:
            pass
        else:
            pass


if __name__ == "__main__":
    pN = PreNetwork(species, reactions)
    ckappas = np.array([0.2, 0.3, 0.4, 0.5]) #in order is I, E, F, C
    assert len(species) == len(X), "species set not aligned with initial population"
    simulate_compartments(pN, sample_pmf_X, ctime_steps, ckappas)


