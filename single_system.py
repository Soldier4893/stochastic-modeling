from reaction_network import PreNetwork, SingleReactionNetwork
import numpy as np


def main(reactions, species, X, time_steps):
    # temp
    np.random.seed(19)
    
    assert len(species) == len(X), "species set not aligned with initial population"
    pN = PreNetwork(species, reactions)
    zeta_table, rate_table, kappas = pN.get_tables()

    network = SingleReactionNetwork(X, zeta_table, rate_table, kappas, time_steps)
    network.simNextReaction()
    print()
    network.simGillespie()
    #network.graph()

if __name__ == "__main__":
    reactions = [
        ("G","G+M", 200),
        ("M","M+P", 10),
        ("M","0", 25),
        ("P","0", 1),
    ]
    # reactions = [
    #     ("G","G+M", 200),
    #     ("M","M+P", 10),
    #     ("M","0", 25),
    #     ("P","0", 1),
    #     ("2P","D", 0.01),
    #     ("D","0", 1),
    # ]
    species = ['G','M','P', 'D']
    X = np.array([1, 0, 0, 0], dtype = np.int32)
    time_steps = 10

    main(reactions, species, X, time_steps)


