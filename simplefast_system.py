from compartment_manager import SimpleFastCompartmentManager, SimpleFastCompartmentManager2, SimpleFastCompartmentManager1
from multiprocessing import Pool
from matplotlib import pyplot as plt
import numpy as np
import time
import cProfile

timelim = 2

def create_simulation():
    x = sfcm.simComparts(timelim)
    # print()
    return np.sum(x)

def create_simulation1():
    y, o = sfcm1.simComparts(timelim)
    # print(z)
    return o
    # return np.sum(y), 

def create_simulation2():
    z, t, o = sfcm2.simComparts(timelim)
    # print('help')
    return o
    # return np.sum(z) + t


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
    # cProfile.run('create_simulation()')
    sum = 0
    sum1 = 0
    sum2 = 0
    timer = 0
    timer1 = 0
    timer2 = 0

    count = 1000
    sfcm2 = SimpleFastCompartmentManager2((1, 1, 2.01, 0), 1, 1)
    sfcm2.simComparts(100)
    sfcm2.graph()

    sfcm1 = SimpleFastCompartmentManager1((1, 1, 2.01, 0), 1, 1)
    sfcm1.simComparts(100)
    sfcm1.graph()

    plt.legend(['1 numS','1 numC', '2 numS','2 numC'])
    plt.show()
    # for _ in range(count):
    #     sfcm2 = SimpleFastCompartmentManager2((1, 1, 2, 0), 1, 1)
    #     sfcm = SimpleFastCompartmentManager((1, 1, 2, 0), 1, 1)
    #     sfcm1 = SimpleFastCompartmentManager1((1, 1, 2, 0), 1, 1)

    #     s = time.time()
    #     sum += create_simulation()
    #     timer +=time.time()-s

    #     s = time.time()
    #     sum1 += create_simulation1()
    #     timer1 +=time.time()-s
        
    #     s = time.time()
    #     sum2+=create_simulation2()
    #     timer2 +=time.time()-s

    # print(timer/count, timer1/count, timer2/count)
    # print(sum, sum1, sum2)
















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

