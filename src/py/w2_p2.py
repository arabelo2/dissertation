# Import libraries
import numpy as np
from math import *
import matplotlib.pyplot as plt

# Initiate parameters
xmax = 10.0                     # physical domain (m)
nx = 100                        # number of space sample
a = 0.25                        # exponent of Gaussian function
dx = xmax/(nx - 1)              # grid spacing dx (m)
x0 = xmax/2                     # center of Gaussian function x0 (m)
x = np.linspace(0, xmax, nx)    # defining space variable

# Initialization of Gaussian function
f = (1./np.sqrt(2*np.py*a))*np.exp(-(((x - x0)**2)/(2*a)))

# Plotting of Gaussian
plt.plot(x, f)
plt.title('Gaussian function')
plt.xlabel('x, m')
plt.ylabel('Amplitude')
plt.xlim(0, xmax)
plt.grid()
plt.show()