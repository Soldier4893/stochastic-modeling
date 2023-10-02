import reaction_network as rn
import numpy as np

def main(reactions, species, X, time_steps):
    assert len(species) == len(X), "species set not aligned with initial population"
    pN = rn.preNetwork(species, reactions)
    zeta_table, rate_table, rates = pN.get_tables()

    network = rn.reactionNetwork(X, zeta_table, rate_table, rates, time_steps)
    network.simulate()
    network.graph()

if __name__ == "__main__":
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
    time_steps = 10000

    main(reactions, species, X, time_steps)


