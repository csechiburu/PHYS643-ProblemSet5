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

# Step 1: create a function that, given 2 arrays f1 and f2, spits out the fluxes J1 and J2
 
# Step 2: set up number of grid points (i.e., number of cell centers), number of time steps, and their sizes

Ngrid = 100 
Ntime = 1000
dx = 1.
dt = 0.01
# Value of speed of sound (constant!)
cs = 1.

# Step 3: set up spatial array x, gaussian perturbation, f1, f2, and corresponding J1, J2

x = np.arange(Ngrid)
# Gaussian with amplitude "a", width "sigma"
a = 10
sigma = 3.0
mu = Ngrid/2.0
gaussian =  a*np.exp(-(x-mu)**2/(2.*sigma**2))

f1 = 1.0 + np.copy(gaussian)
f2 = np.copy(gaussian)

# Step : set up plot (from Eve's code)

plt.ion()
fig, axes = plt.subplots(1,1)
axes.set_title('f1 - density')
# Plot solutions for each time step
plot_f1,= axes.plot(x,f1,'k-')
# Fix axes limits for better visualization
axes.set_xlim([0,Ngrid])
#axes.set_ylim([-0.0,3.0])
fig.canvas.draw()

# Step : solve hydro equations over time

count = 0

while count < Ntime:
    J1 = np.zeros(Ngrid+1)
    J2 = np.zeros(Ngrid+1)
    for i in range(1,Ngrid):
        # Get u array
        u = 0.5*(f2[i-1]/f1[i-1] + f2[i]/f1[i])
        # Compute J1 and J2
        if u > 0:
            J1[i] = u*f1[i-1]
            J2[i] = u*f2[i-1]
        if u < 0:
            J1[i] = u*f1[i]
            J2[i] = u*f2[i]
    # Continuity equation - f1
    f1 = f1 - (dt/dx)*(J1[1:]-J1[:-1])
    # Boundary conditions 
    f1[0] = f1[0]-(dt/dx)*J1[0]
    f1[-1] = f1[-1]+(dt/dx)*J1[-2]
    # Momentum equation - f2
    f2 = f2 - (dt/dx)*(J2[1:]-J2[:-1])
    f2[1:-1] = f2[1:-1] - (dt*(cs**2)/dx)*(f1[2:]-f1[:-2])
    # Boundary conditions 
    f2[0] = f2[0]-(dt/dx)*J2[0]
    f2[-1] = f2[-1]+(dt/dx)*J2[-2]
    # Plotting each
    plot_f1.set_ydata(f1)
    #plot_f2.set_ydata(f2)
    fig.canvas.draw()
    plt.pause(0.001)
    count += dt

''' These are comments for myself :) please ignore
1) continuity/momentum <- need to solve f <- need to solve j <- need to solve u <- u depends on f

-/ define arrays f1 and f2
-/ take f1 and f2, each at j and j+1 => f, f[1:]
-/ compute u at j+1/2 => create a function that spits out J's, given f's
-/ compute u at j-1/2
-/ determine whether u at j+1/2 is > 0 or < 0
-/ determine whether u at j-1/2 is > 0 or < 0
-/ compute J1 at j+1/2
-/ compute j2 at j-1/2
- use those J's, dt, dx and f1 at j to compute next f1 at j (as we've always done)
- repeat all the above for f2

2) add source term -> need to solve f

- in step 1 we get f at j but for a time n+1/2, not n+1
- For continuity mass is conserved so fj at n+1 = fj at n+1/2
- for momentum fj at n+1 = fj at n+1/2 - (dt/dx)*cs^2*( f1,j+1^(n+1/2) - f1,j-1^(n+1/2))
'''

