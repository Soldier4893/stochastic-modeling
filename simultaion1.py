import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

rates = np.array([200, 10, 25, 1])
G = [1]
M = [0]
P = [0]
X = np.array([1, 0, 0])
zetas = [np.array([0, 1, 0]),np.array([0, 0, 1]),np.array([0, -1, 0]),np.array([0, 0, -1])]
time = [0]

for t in range(5000):
    rate = np.dot(rates, np.array([X[0], X[1], X[1], X[2]]))
    delay = np.random.exponential(1.0/rate)
    ps = np.array([X[0] * rates[0],X[1] * rates[1],X[1] * rates[2],X[2] * rates[3]])
    choose = np.random.choice([0, 1, 2, 3], p = ps/np.sum(ps))

    X += zetas[choose]

    G.append(X[0])
    M.append(X[1])
    P.append(X[2])
    time.append(time[-1] + delay)

plt.plot(time, G)
plt.plot(time, M)
plt.plot(time, P)
plt.legend(['G', 'M', 'P'])
plt.show()


