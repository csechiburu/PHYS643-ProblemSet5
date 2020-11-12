'''
Problem #1
----------
Solving Advection equation using FTCS and Lax-Friedrichs methods
----------
@author: Constanza Echiburu
Nov. 12th 2020
'''

# imports
import numpy as np
import matplotlib.pyplot as plt


# Step 1: set up number of grid points, number of time steps, and sizes of steps

Ngrid = 50
Ntime = 1000
dx = 1
dt = 1

# Step 2: define known variables

u = -0.1
alpha = u*dt/(2.*dx)

# Step 3: check if chosen variable meet Courant condition:
print('Courant condition is',dt<=abs(dx/u))

# Step 4: set up spatial array over which f is defined, with initial conditions

x = np.arange(Ngrid) 
# Array for FTCS
f = np.copy(x)*(1./Ngrid)
# Array for Lax-Friedrich
f1 = np.copy(x)*(1./Ngrid)


# Step 5: set up plot (from Eve's code)

plt.ion()
fig, axes = plt.subplots(1,2)
axes[0].set_title('FTCS')
axes[1].set_title('Lax-Friedrichs')
# Plot initial state for reference
axes[0].plot(x,f,'k-')
axes[1].plot(x,f1,'k-')
# Plot solutions for each time step
plot_ftcs,= axes[0].plot(x,f,'ro')
plot_lf,= axes[1].plot(x,f1,'ro')
# Fix axes limits for better visualization
axes[0].set_xlim([0,Ngrid])
axes[0].set_ylim([-0.1,2.0])
axes[0].set_xlabel('x')
axes[0].set_ylabel('f')
axes[1].set_xlim([0,Ngrid])
axes[1].set_ylim([-0.1,2.0])
axes[1].set_xlabel('x')
axes[1].set_ylabel('f')
fig.canvas.draw()

# Step 6: solve Advection equation over time

count = 0

while count < Ntime:
    # FTCS
    f[1:-1] = f[1:-1] - alpha*(f[2:] - f[:-2])
    # Lax-Friedrich
    f1[1:-1] = 0.5*(f1[2:] + f1[:-2]) - alpha*(f1[2:] - f1[:-2])
    # Plotting each
    plot_ftcs.set_ydata(f)
    plot_lf.set_ydata(f1)
    fig.canvas.draw()
    plt.pause(0.001)
    count += dt






