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
f = (1./np.sqrt(2*np.pi*a))*np.exp(-(((x - x0)**2)/(2*a)))

# Plotting of Gaussian
plt.plot(x, f)
plt.title('Gaussian function')
plt.xlabel('x, m')
plt.ylabel('Amplitude')
plt.xlim(0, xmax)
plt.grid()
plt.show()

# Second derivative with 3-point operator

# Initiation of numerical and analytical derivatives
nder3 = np.zeros(nx)            # numerical derivative
ader = np.zeros(nx)             # analytical derivative

# Numerical second derivative of the given function
for i in range(1, nx - 1):
    nder3[i] = (f[i + 1] - 2*f[i] + f[i - 1])/(dx**2)

# Analytical second derivative of the Gaussian function
ader = 1./np.sqrt(2*np.pi*a)*((x - x0)**2/a**2 -1/a)*np.exp(-1/(2*a)*(x - x0)**2)

# Exclude boundaries
ader[0] = 0.
ader[nx - 1] = 0.

# Calculate rms error of numerical derivative
rms = np.sqrt(np.mean((nder3 - ader)**2))

# Plotting
plt.plot(x, nder3, label="Numerical derivative, 3 points", lw=2, color="violet")
plt.plot(x, ader, label="Analytical derivative", lw=2, ls="--")
plt.plot(x, nder3 - ader, label="Difference", lw=2, ls=":")
plt.title("Second derivative, Err (rms) = %.6f" % (rms))
plt.xlabel('x, m')
plt.ylabel('Amplitude')
plt.legend(loc="lower left")
plt.grid()
plt.show()

# First derivative with 5-point operator

# Initiation of numerical derivative
nder5 = np.zeros(nx)            # numerical derivative

# Calculation of second derivative of the given function
for i in range(1, nx - 2):
    nder5[i] = (-1./12*f[i - 2] + 4./3*f[i - 1] - 5./2*f[i] + 4./3*f[i + 1] - 1./12*f[i + 2])/(dx**2)

# Exclude boundaries
ader[1] = 0.
ader[nx - 2] = 0.

# Calculate rms error of numerical derivative
rms = rms*0
rms = np.sqrt(np.mean((nder5 - ader)**2))

# Plotting
plt.plot(x, nder5, label="Numerical derivative, 5 points", lw=2, color="violet")
plt.plot(x, ader, label="Analytical derivative", lw=2, ls="--")
plt.plot(x, nder5 - ader, label="Difference", lw=2, ls=":")
plt.title("Second derivative, Err (rms) = %.6f" % (rms))
plt.xlabel('x, m')
plt.ylabel('Amplitude')
plt.legend(loc="lower left")
plt.grid()
plt.show()