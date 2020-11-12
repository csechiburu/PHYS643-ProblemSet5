'''
Problem #2
----------
Solving Diffusion & Advection-Difussion equations
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

# Step 2: define bulk velocity, alpha, diffusion coefficients and the corresponding betas

u = -0.1
alpha = u*dt/dx
D1, D2 = 0.1, 1.0
beta1, beta2 = D1*dt/dx**2, D2*dt/dx**2


# Step 3: check if chosen variable meet Courant condition:
print('Courant condition is',dt<=abs(dx/u))

# Step 4: set up spatial array over which f is defined, with initial conditions

x = np.arange(Ngrid) 
f1 = np.copy(x)*(1./Ngrid)
f2 = np.copy(x)*(1./Ngrid)

# Step 5: build matrix A, where first term is the diagonal, and the other two are off diagonal

A1 = np.eye(Ngrid)*(1.0+2.0*beta1) + np.eye(Ngrid,k=-1)*(-beta1) + np.eye(Ngrid,k=1)*(-beta1)
A2 = np.eye(Ngrid)*(1.0+2.0*beta2) + np.eye(Ngrid,k=-1)*(-beta2) + np.eye(Ngrid,k=1)*(-beta2)

# Ensure boundary conditions
A1[0][0] = 1.
A1[0][1] = 0
A1[-1][-1] = 1. + beta1
A2[0][0] = 1.
A2[0][1] = 0
A2[-1][-1] = 1. + beta2

# Step 6: set up plot (from Eve's code)

plt.ion()
fig, axes = plt.subplots(1,2)
axes[0].set_title('D='+str(D1))
axes[1].set_title('D='+str(D2))
# Plot initial state for reference
axes[0].plot(x,f1,'k-')
axes[1].plot(x,f2,'k-')
# Plot solutions for each time step
plot_D1,= axes[0].plot(x,f1,'ro')
plot_D2,= axes[1].plot(x,f2,'ro')
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

# Step 7: solve Diffusion equation over time

count = 0

while count < Ntime:
    # Diffusion term
    f1[1:-1] = np.linalg.solve(A1,f1)[1:-1]
    f2[1:-1] = np.linalg.solve(A2,f2)[1:-1]
    # Advection terms
    f1[1:-1] = 0.5*(f1[2:]+f1[:-2]) - 0.5*alpha*(f1[2:]-f1[:-2])
    f2[1:-1] = 0.5*(f2[2:]+f2[:-2]) - 0.5*alpha*(f2[2:]-f2[:-2])
    # Plotting each
    plot_D1.set_ydata(f1)
    plot_D2.set_ydata(f2)
    fig.canvas.draw()
    plt.pause(0.001)
    count += dt

