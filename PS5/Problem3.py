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

def get_J(f1,f2):
    ''' A function to get the fluxes at the interface between the donor cells.
    Input: f1 - array, to represent density
           f2 - array, to represent density*velocity
    Output: J1_r, J1_l - arrays, to represent flux of f1 at right/left interface
            J2_r, J2_l - arrays, to represent flux of f2 at right/left interface
    '''
    # Get u array at j+1/2 (to the right)
    u_r = 0.5*(f2[1:-1]/f1[1:-1]+f2[2:]/f1[2:])
    # Get u array at j-1/2 (to the left)
    u_l = 0.5*(f2[1:-1]/f1[1:-1]+f2[:-2]/f1[:-2])
    # Create J1 and J2 arrays, same length as u
    J1_r = np.arange(len(u_r))
    J2_r = np.arange(len(u_r))
    J1_l = np.arange(len(u_r))
    J2_l = np.arange(len(u_r))
    # Checking signs of u
    for j in range(len(J1_r)):
        # Compute J1 and J2 at j+1/2
        if u_r[j] > 0:
            J1_r[j] = u_r[j]*f1[j]
            J2_r[j] = u_r[j]*f2[j]
        if u_r[j] < 0:
            J1_r[j] = u_r[j]*f1[j+1]
            J2_r[j] = u_r[j]*f2[j+1]
        # Compute J1 and J2 at j-1/2
        if u_l[j] > 0:
            J1_l[j] = u_l[j]*f1[j-1]
            J2_l[j] = u_l[j]*f2[j-1]
        if u_l[j] < 0:
            J1_l[j] = u_l[j]*f1[j]
            J2_l[j] = u_l[j]*f2[j]       
    return J1_r, J1_l, J2_r, J2_l
    
# Step 2: set up number of grid points (i.e., number of cell centers), number of time steps, and their sizes

Ngrid = 50 
Ntime = 5000
dx = 1
dt = 1
# Value of speed of sound (constant!)
cs = 0.1

# Step 3: set up spatial array x, gaussian perturbation, f1, f2, and corresponding J1, J2

x = np.arange(Ngrid)
# Gaussian with amplitude "a", width "sigma"
a = 1
sigma = 4.0
gaussian =  a*np.exp(-(x-25)**2/(2*sigma**2))

f1 = 1.0 + np.copy(gaussian)
f2 = np.copy(gaussian)
J1_r, J1_l, J2_r, J2_l = get_J(f1,f2)

# Step 4: Boundary conditions 
f1[0] = f1[0]-(dt/dx)*J1_l[0]
f1[-1] = f1[-1]+(dt/dx)*J1_r[-2]
f2[0] = f2[0]-(dt/dx)*J2_l[0]
f2[-1] = f2[-1]+(dt/dx)*J2_r[-2]

# Step : set up plot (from Eve's code)

plt.ion()
fig, axes = plt.subplots(1,2)
axes[0].set_title('f1')
axes[1].set_title('f2')
# Plot solutions for each time step
plot_f1,= axes[0].plot(x,f1,'-')
plot_f2,= axes[1].plot(x,f2,'-')
# Fix axes limits for better visualization
axes[0].set_xlim([0,Ngrid])
axes[0].set_ylim([-5.0,5.0])
axes[1].set_xlim([0,Ngrid])
axes[1].set_ylim([-5.0,5.0])
fig.canvas.draw()

# Step : solve hydro equations over time

count = 0

while count < Ntime:
    # Continuity equation - f1
    f1[1:-1] = f1[1:-1] - (dt/dx)*(J1_r-J1_l) 
    # Momentum equation - f2
    f2[1:-1] = (f2[1:-1] - (dt/dx)*(J2_r-J2_l)) - (dt*(cs**2)/dx)*(f1[2:]-f1[:-2])  
    # Plotting each
    plot_f1.set_ydata(f1)
    plot_f2.set_ydata(f2)
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

