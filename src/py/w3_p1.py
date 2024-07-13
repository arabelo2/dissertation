# Import Libraries
import numpy as np
import matplotlib.pyplot as plt

# Sub-plot Configuration
# ----------------------
from matplotlib import gridspec

# Ignore Warning Messages
# -----------------------
import warnings
warnings.filterwarnings("ignore")

# Paramenter configuration
nx = 10000                      # number of grid points in x-direction
xmax = 10000                    # physical domain (m)
dx = xmax/(nx - 1)              # grid point distance in x-direction
c0 = 334.                       # wave speed in medium (m/s)
isrc = int(nx/2)                # source location in grid in x-direction
nt = 1000                       # maximum number of time steps
dt = 0.002                     # time step

# Source time function parameters
f0 = 25                         # dominant frequency of the source (Hz)
t0 = 4./f0                      # source time shift

# Snapshot
idisp = 5                       # display frequency

# Source time function (Gaussian)
src = np.zeros(nt + 1)
time = np.linspace(0*dt, nt*dt, nt)

# First derivative of a Gaussian
# W3V5 Python: Waves in 1D (Comment or uncomment)
src = -2.*(time - t0)*(f0**2)*(np.exp(-1.*(f0)**2*(time - t0)**2))
# W3P1 Acoustic Waves 1D (Lab) (Comment or uncomment)
# src = -8.*(time - t0)*(f0**2)*(np.exp(-1.*(4*f0)**2*(time - t0)**2))

# Plot source time function

# Plot position configuration

with plt.ion():
    fig1 = plt.figure(figsize=(10, 6))
    gs1 = gridspec.GridSpec(1, 2, width_ratios=[1, 1], hspace=0.3, wspace=0.3)

    # Plot source time function
    ax1 = plt.subplot(gs1[0])
    ax1.plot(time, src)             # plot source time function
    ax1.set_title('Source time function')
    ax1.set_xlim(time[0], time[-1])
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Amplitude')

    # Plot source spectrum
    ax2 = plt.subplot(gs1[1])
    spec = np.fft.fft(src)  # source time function in frequency domain
    freq = np.fft.fftfreq(spec.size, d=dt)  # time domain to frequency domain
    ax2.plot(np.abs(freq), np.abs(spec))  # plot frequency and amplitude
    ax2.set_xlim(0, 250)  # only display frequency from 0 to 250 Hz
    ax2.set_title('Source spectrum')
    ax2.set_xlabel('Frequency (Hz)')
    ax2.set_ylabel('Amplitude')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    plt.show()

# Plot snapshot & seismogran

# Initialize empty pressure
p = np.zeros(nx)            # p at time n (now)
pold = np.zeros(nx)         # p at time n - 1 (past)
pnew = np.zeros(nx)         # p at time n + 1 (present)
d2px = np.zeros(nx)         # 2nd space derivative of p

# Initialize model (assume homogeneous model)
c = np.zeros(nx)
c = c + c0                  # initialize wave velocity in model

# Initialize coordinate
x = np.arange(nx)
x = x*dx                  # coordinate in x-direction

# plot position configuration
with plt.ion():
    fig2 = plt.figure(figsize=(10, 6))
    gs2 = gridspec.GridSpec(1, 1, width_ratios=[1], hspace=0.3, wspace=0.3)

    # Plot 1D wave propagation
    # Note: comma is needed to update the variable
    ax3 = plt.subplot(gs2[0])

    # Plot position of the source in snapshot
    leg1, = ax3.plot(isrc, 0, 'r*', markersize=11)

    # Plot position of the receiver in snapshot
    # leg2, = ax3.plot(ir, 0, 'k^', markerzise=8)

    # Plot pressure update each time step
    up31, = ax3.plot(p)

    ax3.set_xlim(0, xmax)
    ax3.set_ylim(-np.max(p), np.max(p))
    ax3.set_title('Time step (nt) = 0')
    ax3.set_xlabel('x (m)')
    ax3.set_ylabel('Pressure Amplitude')
    # ax3.legend((leg1, leg2), ('Source', 'Receiver'), loc='upper right',\
    #     fontsize=10, numpoints=1)
    plt.show()

# 1D Wave Propagation (Finite Difference Solution)

# Loop over time
for it in range(nt):

    # 2nd derivative in space
    for i in range(1, nx - 1):
        d2px[i] = (p[i + 1] - 2 * p[i] + p[i - 1]) / dx ** 2

    # Time Extrapolation
    pnew = 2 * p - pold + c ** 2 * dt ** 2 * d2px

    # Add Source Term at isrc
    # Absolute pressure w.r.t analytical solution
    pnew[isrc] = pnew[isrc] + src[it] / (dx) * dt ** 2

    # Remap Time Levels
    pold, p = p, pnew

    # Plot pressure field
    if (it % idisp) == 0:
        ax3.set_title('Time Step (nt) = %d' % it)
        ax3.set_ylim(-1.1*np.max(abs(p)), 1.1*np.max(abs(p)))

        # plot around propagating wave
        window = 100
        xshift = 25
        ax3.set_xlim(isrc*dx + c0*it*dt - window*dx - xshift, isrc*dx + c0*it*dt + window*dx - xshift)
        up31.set_ydata(p)
        plt.gcf().canvas.draw()
        plt.pause(.001)
