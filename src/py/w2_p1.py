# Import libraries
import numpy as np
from math import *
import matplotlib.pyplot as plt

# Initiate parameters
xmax = 10.0                     # physical domain (m)
nx = 200                        # number of space sample
dx = xmax/(nx - 1)              # grid spacing dx (m)
x = np.linspace(0, xmax, nx)    # defining space variable

# Initialization of sin function
l = 20*dx                       # wavelength
k = 2*np.pi/l                   # wavenumber
f = np.sin(k*x)

# Plot sin function
plt.plot(x, f)
plt.title('Sin function')
plt.xlabel('x, m')
plt.ylabel('Amplitude')
plt.xlim(0, xmax)
plt.grid()
plt.show()

# First derivative with two points

# Initiation of numerical and analytical derivatives
nder = np.zeros(nx)             # numerical derivative
ader = np.zeros(nx)             # analytical derivative

# Numerical derivative of the given function
for i in range(1, nx - 1):
    nder[i] = (f[i + 1] - f[i - 1])/(2*dx)

# Analytical derivative of the given function
ader = k*np.cos(k*x)

# Exclude boundries
ader[0] = 0.
ader[nx -1] = 0.

# Error (rms)
rms = np.sqrt(np.mean(nder - ader)**2)

# Plotting
plt.plot(x, nder, label="Numerical derivative, 2 points", marker='+', color="blue")
plt.plot(x, ader, label="Analytical derivative", lw=2, ls="-", color="black")
plt.plot(x, nder-ader, label="Difference", lw=2, ls=":")
plt.title("First derivative, Err (rms) = %.6f " % (rms))
plt.xlabel('x, m')
plt.ylabel('Amplitude')
plt.legend(loc='lower left')
plt.grid()
plt.show()

# Plotting number of points per wavelength

plt.plot(x, nder, label="Number derivative, 2 points", marker='+', color="blue")
plt.title("First derivative, Error = %.6f, $n_\\lambda$ = %.2f" % (rms, l/dx))
plt.xlabel('x, m')
plt.ylabel('Amplitude')
plt.legend(loc='lower left')
plt.xlim(xmax/2 - 1, xmax/2 + 1)
plt.grid()
plt.show()

# Define the range of number of points per wavelength: [nmin = 4, 5, 6, ..., nmax = 40].
# Loop over points, calculate corresponding wavelength and calculate error.

# Initialize vectors
nmin = 3
nmax = 16
na = np.zeros(nmax - nmin + 1)              # vector with number of points per wavelength
err = np.zeros(nmax - nmin + 1)            # vector with error
j = -1                                      # array index

# Loop through finite-difference derivative calculation
for n in range(nmin, nmax):
    j = j + 1                               # array index
    na[j] = n
    
    # Initialize sin function
    l = na[j]*dx                            # wavelength
    k = 2*np.pi/l                           # wavenumber
    f = np.sin(k*x)
    
    # Numerical derivative of the sin function
    for i in range(1, nx - 1):
        nder[i] = (f[i + 1] - f[i - 1])/(2*dx)
        
    # Analytical derivative of the sin function
    ader = k*np.cos(k*x)
    
    # Exclude boundaries
    ader[0] = 0.
    ader[nx - 1] = 0.
    
    i0 = np.int64(nx/2)
    # Error (rms)
    err[j] = ((nder[i0] - ader[i0])**2/ader[i0]**2)*100
    
# Plotting error as function of number of points per wavelength
plt.plot(na, err, ls='-', color="blue")
plt.title('Error as a function of n$_\\lambda$')
plt.xlabel('n$_\\lambda$')
plt.ylabel('rms')
plt.grid()
plt.show()