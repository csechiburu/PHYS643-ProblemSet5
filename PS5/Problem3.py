'''
Problem #3
----------
1-D Hydro solver
----------
@author: Constanza Echiburu
Nov. 12th 2020
'''

# imports
import numpy as np
import matplotlib.pyplot as plt
 
# Step 1: set up number of grid points (i.e., number of cell centers), number of time steps, and their sizes

Ngrid = 100 
Ntime = 1000
dx = 1.
dt = 0.08
# Value of speed of sound (constant!)
cs = 1.

# Step 2: set up spatial array x, gaussian perturbation, f1, and f2

x = np.arange(Ngrid)
# Gaussian with amplitude "a", width "sigma", centered at "mu"
a = 1
sigma = 3.0
mu = Ngrid/2.0
gaussian =  a*np.exp(-(x-mu)**2/(2.*sigma**2))
# density and density*velocity arrays
f1 = 1.0 + np.copy(gaussian)
f2 = np.copy(gaussian)

# Step 3: set up plot (from Eve's code)

plt.ion()
fig, axes = plt.subplots(1,1)
axes.set_title('Density Evolution')
axes.set_ylabel('f1')
axes.set_xlabel('x')
# Plot solutions for each time step
plot_f1,= axes.plot(x,f1,'k-')
# Fix axes limits for better visualization
axes.set_xlim([0,Ngrid])
#axes.set_ylim([-0.0,3.0])
fig.canvas.draw()

# Step 4: solve hydro equations over time

count = 0

while count < Ntime:
    # Create J arrays to store the flux values
    J1 = np.zeros(Ngrid+1)
    J2 = np.zeros(Ngrid+1)
    for i in range(1,Ngrid):
        # Get u array
        u = 0.5*(f2[i-1]/f1[i-1] + f2[i]/f1[i])
        # Compute J1 and J2 according to the sign of u
        if u > 0:
            J1[i] = u*f1[i-1]
            J2[i] = u*f2[i-1]
        if u < 0:
            J1[i] = u*f1[i]
            J2[i] = u*f2[i]
    # Compute continuity equation - f1
    f1 = f1 - (dt/dx)*(J1[1:]-J1[:-1])
    # Boundary conditions 
    f1[0] = f1[0]-(dt/dx)*J1[0]
    f1[-1] = f1[-1]+(dt/dx)*J1[-2]
    # Compute Euler equation without source term - f2
    f2 = f2 - (dt/dx)*(J2[1:]-J2[:-1])
    # now add source term (pressure gradient)
    f2[1:-1] = f2[1:-1] - (dt*(cs**2)/dx)*(f1[2:]-f1[:-2])
    # Boundary conditions 
    f2[0] = f2[0]-(dt/dx)*J2[0]
    f2[-1] = f2[-1]+(dt/dx)*J2[-2]
    # Plotting
    plot_f1.set_ydata(f1)
    fig.canvas.draw()
    plt.pause(0.001)
    count += dt

