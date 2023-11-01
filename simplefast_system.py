from compartment_manager import SimpleFastCompartmentManager, SimpleFastCompartmentManager2, SimpleFastCompartmentManager1
from multiprocessing import Pool
from matplotlib import pyplot as plt
import numpy as np
import time
import cProfile

timelim = 20
# np.random.seed(19)

# def create_simulation():
#     x = sfcm.simComparts(timelim)
#     # print()
#     return np.sum(x)

def create_simulation1():
    c, s = sfcm1.simComparts(timelim)
    # print(z)
    return c, s

def create_simulation2():
    c, s = sfcm2.simComparts(timelim)
    # print('help')
    return c, s


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
    sums = [0,0,0,0]
    timer1 = 0
    timer2 = 0

    # count = 10
    sfcm1 = SimpleFastCompartmentManager1((1, 1, 2, 0), 1, 1)
    sfcm1.simComparts(2500)
    sfcm1.graph()
    
    # sfcm2 = SimpleFastCompartmentManager2((1, 1, 2.5, 0), 1, 1)
    # sfcm2.simComparts(5)
    # sfcm2.graph()

    # plt.legend(['1 numS','1 numC','2 numS','2 numC'])
    # plt.show()
    count = 50
    # for _ in range(count):
    #     # sfcm2 = SimpleFastCompartmentManager2((1, 1, 1.5, 0), 1, 1)
    #     # sfcm = SimpleFastCompartmentManager((1, 1, 2, 0), 1, 1)
    #     sfcm1 = SimpleFastCompartmentManager1((1, 1, 2, 0), 1, 1)

    #     s = time.time()
    #     a, b = create_simulation1()
    #     sums[0] += a
    #     sums[1] += b
    #     timer1 +=time.time()-s
        
    #     # s = time.time()
    #     # a, b = create_simulation2()
    #     # sums[2] += a
    #     # sums[3] += b
    #     # timer2 +=time.time()-s

    # print(timer1/count, timer2/count)
    # print("1c", sums[0]/count)
    # print("1s", sums[1]/count)
    # print("2c", sums[2]/count)
    # print("2s", sums[3]/count)
















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

