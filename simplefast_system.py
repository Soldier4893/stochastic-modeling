from compartment_manager import SimpleFastCompartmentManager
from multiprocessing import Pool
import numpy as np
import time

def create_simulation(t):
    sfcm = SimpleFastCompartmentManager((1, 1, 3.1, 0), 1, 1)
    x = sfcm.simComparts(t)
    print(x)

def test(threads, duration, many):
    s = time.time()
    with Pool(threads) as p:
        p.map(create_simulation, [duration for i in range(many)])
    print(threads,': ', duration,' for ',many,' ',time.time()-s)

long = 120
short = 30
many = 150
few = 40
threads = 8
threads2 = 4

if __name__ == "__main__":
    create_simulation(20)
    print('done')
    # test(threads, short, many)
    # test(threads2, short, many)
    # test(threads, long ,many)
    # test(threads2, long ,many)
    # print()
    # test(threads, long ,few)
    # test(threads2, long ,few)
    # print()
    # test(threads, short ,many)
    # test(threads2, short ,few)
    # print()
    # test(threads2, short ,many)
    # test(threads, short ,few)

