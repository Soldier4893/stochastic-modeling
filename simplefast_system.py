from compartment_manager import SimpleFastCompartmentManager
from multiprocessing import Pool
import numpy as np
import time

def create_simulation(t):
    sfcm = SimpleFastCompartmentManager((2, 0.1, 1.3, 0.1), 3, 0.2)
    x = sfcm.simComparts(t)

if __name__ == "__main__":
    s = time.time()
    with Pool(16) as p:
        p.map(create_simulation, [100 for i in range(256)])
    e = time.time()
    print(e-s)

    s = time.time()
    create_simulation(100)
    e = time.time()
    print(e-s)


