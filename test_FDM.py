import math
import numpy as np
from matplotlib import pyplot as plt

count = 10

def f(x):
    return math.sin(x)

def u(h):
    # return (-25*f(0*h)+48*f(1*h)-36*f(2*h)+16*f(3*h)-3*f(4*h))/(12*h)
    
    return (f(-2*h)-8*f(-h)+8*f(h)-f(2*h))/(12*h)
    
def v(h):
    # return (f(h)-f(-h))/(2*h)
    # return (1000*f(-1.1*h)-1331*f(-h)+1331*f(h)-1000*f(1.1*h))/(462*h)
    # return (-f(-3*h)+9*f(-2*h)-45*f(-1*h)+45*f(1*h)-9*f(2*h)+1*f(3*h))/(60*h)
    return (f(-3*h)-27*f(-1*h)+27*f(h)-1*f(3*h))/(48*h)
    # return (f(-3*h)-f(-1*h)+f(h)-f(3*h))/(48*h)
    # return (-11*f(0)+18*f(h)-9*f(2*h)+2*f(3*h))/(6*h)
    # return (-11*f(0)+18*f(h)-9*f(2*h)+2*f(3*h))/(6*h)
    pass

x = (0.5)**np.arange(count)
A = np.ones((2, count))
A[1] = np.log(x)
A = A.T

b = np.log(np.abs(np.vectorize(u)(x) - np.ones(count)))
Y = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, b))
plt.plot(np.log(x), b)
print(Y)


x = (0.5)**np.arange(count)
A = np.ones((2, count))
A[1] = np.log(x)
A = A.T

b = np.log(np.abs(np.vectorize(v)(x) - np.ones(count)))
Y = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, b))
plt.plot(np.log(x), b)
print(Y)

plt.legend(['u, standard','v, experimental'])
plt.show()





