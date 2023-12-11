from compartment_model_sensitivity import SimpleFastCompartmentManager, CoupledCompartments
from multiprocessing import Pool
from matplotlib import pyplot as plt
import numpy as np

def compute_variance_over_time(x1, x2):
    count1 = 0
    count2 = 0
    l1, l2 = len(x1), len(x2)
    data = np.zeros((l1 + l2, 2))

    for i in range(l1 + l2):
        index = count1 + count2
        if count1 == l1:
            data[index, 0] = x2[count2, 0]
            data[index, 1] = abs(x2[count2, 1] - x1[-1, 1])
            count2 += 1
        elif count2 == l2:
            data[index, 0] = x1[count1, 0]
            data[index, 1] = abs(x2[-1, 1] - x1[count1, 1])
            count1 += 1
        else:
            if x1[count1, 0] < x2[count2, 0]:
                data[index, 0] = x1[count1, 0]
                data[index, 1] = abs(x2[count2-1, 1] - x1[count1, 1])
                count1 += 1
            else:
                data[index, 0] = x2[count2, 0]
                data[index, 1] = abs(x2[count2, 1] - x1[count1-1, 1])
                count2 += 1
    return data.T

def plot_variances(time_limit):
    comp = SimpleFastCompartmentManager((1, 1, 1.5 - h, 1), 1, 1)
    comp.simComparts(time_limit)
    C,R,S,P1 = comp.getStats()

    comp = SimpleFastCompartmentManager((1, 1, 1.5 + h, 1), 1, 1)
    comp.simComparts(time_limit)
    C,R,S,P2 = comp.getStats()
    x = compute_variance_over_time(P1, P2)
    plt.plot(x[0], x[1])
    plt.show()

def mc(h, time_limit, num_iterations):
    samplesl1 = np.zeros(num_iterations)
    samplesr1 = np.zeros(num_iterations)

    for i in range(num_iterations):
        np.random.seed(np.random.randint(0, num_iterations * 10))

        comp = SimpleFastCompartmentManager((1, 1, 1.5 - h, 1), 1, 1)
        comp.simComparts(time_limit)
        samplesl1[i]=comp.numC
        del comp

        comp = SimpleFastCompartmentManager((1, 1, 1.5 + h, 1), 1, 1)
        comp.simComparts(time_limit)
        samplesr1[i]=comp.numC
        del comp
    return (samplesr1-samplesl1)/(2*h*num_iterations)

def crn(h, time_limit, num_iterations):
    samplesl2 = np.zeros(num_iterations)
    samplesr2 = np.zeros(num_iterations)

    np.random.seed()
    for i in range(num_iterations):
        comp = SimpleFastCompartmentManager((1, 1, 1.5 - h, 1), 1, 1)
        comp.simComparts(time_limit)
        samplesl2[i]=comp.numC
        del comp

        comp = SimpleFastCompartmentManager((1, 1, 1.5 + h, 1), 1, 1)
        comp.simComparts(time_limit)
        samplesr2[i]=comp.numC
        del comp
    return  (samplesr2-samplesl2)/(2*h*num_iterations)

def cfd(h, time_limit, num_iterations):
    samplesl1 = np.zeros(num_iterations)
    samplesr1 = np.zeros(num_iterations)
    for i in range(num_iterations):
        coupledCompart = CoupledCompartments((1, 1, 1.5, 1), 1, 1)
        coupledCompart.simComparts(time_limit=time_limit, h=h)
        samplesl1[i] = coupledCompart.numCs[0]
        samplesr1[i] = coupledCompart.numCs[1]
    return (samplesr1-samplesl1)/(2*h*num_iterations)

if __name__ == "__main__":
    num_iterations = 1000
    h=0.01
    time_limit = 5
    y = mc(h = h, time_limit = time_limit, num_iterations=num_iterations)
    print('monte carlo:',np.mean(y), ' pm ',np.std(y),'\n')

    z = mc(h = h, time_limit = time_limit, num_iterations=num_iterations)
    print('common rand:',np.mean(z), ' pm ',np.std(z),'\n')

    x = cfd(h = 0.01, time_limit = time_limit, num_iterations=num_iterations)
    print('coupled:',np.mean(x), ' pm ',np.std(x),'\n')




