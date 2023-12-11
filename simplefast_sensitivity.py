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

def plot_variances(timelength):
    comp = SimpleFastCompartmentManager((1, 1, 1.5 - h, 1), 1, 1)
    comp.simComparts(timelength)
    C,R,S,P1 = comp.getStats()

    comp = SimpleFastCompartmentManager((1, 1, 1.5 + h, 1), 1, 1)
    comp.simComparts(timelength)
    C,R,S,P2 = comp.getStats()
    x = compute_variance_over_time(P1, P2)
    plt.plot(x[0], x[1])
    plt.show()



def mc_vs_crn(h, timelength, num_samples, num_iterations):
    samplesl1 = np.zeros(num_iterations)
    samplesr1 = np.zeros(num_iterations)
    samplesl2 = np.zeros(num_iterations)
    samplesr2 = np.zeros(num_iterations)

    for i in range(num_iterations):
        np.random.seed(np.random.randint(0, 100))

        comp = SimpleFastCompartmentManager((1, 1, 1.5 - h, 1), 1, 1)
        comp.simComparts(timelength)
        samplesl1[i]=comp.numC
        del comp

        comp = SimpleFastCompartmentManager((1, 1, 1.5 + h, 1), 1, 1)
        comp.simComparts(timelength)
        samplesr1[i]=comp.numC
        del comp

    np.random.seed()
    for i in range(num_iterations):
        comp = SimpleFastCompartmentManager((1, 1, 1.5 - h, 1), 1, 1)
        comp.simComparts(timelength)
        samplesl2[i]=comp.numC
        del comp

        comp = SimpleFastCompartmentManager((1, 1, 1.5 + h, 1), 1, 1)
        comp.simComparts(timelength)
        samplesr2[i]=comp.numC
        del comp

    crn = (samplesr1-samplesl1)/(2*h*num_iterations)
    nmc = (samplesr2-samplesl2)/(2*h*num_iterations)
    return np.mean(crn), np.std(crn), np.mean(nmc), np.std(nmc)

def cfd(num_iterations, step_size):
    samplesl1 = np.zeros(num_iterations)
    samplesr1 = np.zeros(num_iterations)
    for i in range(num_iterations):
        coupledCompart = CoupledCompartments((1, 1, 1.5, 1), 1, 1)
        coupledCompart.simComparts(time_limit=5, h=step_size)
        samplesl1[i] = coupledCompart.numCs[0]
        samplesr1[i] = coupledCompart.numCs[1]
    return (samplesr1[i]-samplesl1[i])/(2*step_size*num_iterations)

if __name__ == "__main__":
    # print(mc_vs_crn(h = 0.01, timelength = 5, num_samples = 60, num_iterations=100))

    # coupledCompart = CoupledCompartments((1, 1, 1.5, 1), 1, 1)
    # coupledCompart.simComparts(time_limit=500, h=0.01)
    # print(coupledCompart.numCs[0], coupledCompart.numCs[1])
    x = cfd(1000, 1)
    print(np.mean(x), ' pm ',np.std(x))




