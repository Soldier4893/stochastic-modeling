from reaction_network import PreNetwork, SingleReactionNetwork, ReactionNetwork
import numpy as np
import multiprocessing


def temp(time):
    reactions = [
        ("G","G+M", 200),
        ("M","M+P", 10),
        ("M","0", 25),
        ("P","0", 16),
    ]
    species = ['G','M','P']
    X = np.array([1, 0, 0], dtype = np.int32)
    pN = PreNetwork(species, reactions)
    zeta_table, rate_table, kappas = pN.get_tables()
    network = ReactionNetwork(X, zeta_table, rate_table, kappas)
    return network.simGillespie(time)

def main():
    with multiprocessing.Pool(4) as p:
        results = p.map(temp, [5 for i in range(200)])
    avg = sum(results) / len(results)
    print(avg)

if __name__ == "__main__":
    main()


