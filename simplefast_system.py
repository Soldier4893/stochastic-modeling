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

def numReaction_vs_time_logfit(count, gap):
    outs = np.zeros(count)
    for i in range(1, count+1):
        sfcm1 = SimpleFastCompartmentManager1((1, 1, 2.1, 0), 1, 1)
        x,y,z =sfcm1.simComparts(i*gap)
        outs[i-1] = z
    A = np.ones((2, count))
    A[1] = np.log(gap*np.arange(1, count+1))
    A = A.T
    b = np.log(np.array(outs))
    x = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, b))
    print(x)
    plt.plot(np.log(gap*np.arange(1, count+1)), np.log(np.array(outs)))
    plt.show()
    return x # K, p


if __name__ == "__main__":
    # cProfile.run('create_simulation()')
    # k, p = numReaction_vs_time_logfit(5, 30)
    

    

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

