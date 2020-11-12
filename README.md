# PHYS643-ProblemSet5
------------------------------------
1. Name: Constanza Echiburu
------------------------------------
2. Version of Python: 3.x
------------------------------------
3. List of python code files and what they do:

- Problem1.py

This code solves the Advection equation for a quantity "f", defined over a spatial array "x", and which evolves with time. It uses the initial condition f(x,t=0)=x over all x, and bulk velocity u=-0.1.

The equation is solved using two different methods: FTCS (Forward-Time Central-Space), and Lax-Friedrichs.

The outputs are two plots showing f(x,t) solved with FTCS (left panel), which is numerically unstable, and Lax-Friedrichs (right panel) which is stable. Numerical stability is set by the Courant condition according to spatial and time steps, and the bulk velocity.


- Problem2.py

This code solves the Advection-Diffusion equation, using the Lax-Friedrich method for advection and the implicit method for diffusion.

It uses the same conditions as in Problem 1, and the output are two plots showing f(x,t). The left panel displays the evolution for constant of diffusivity D=0.1, while the right panel displays the evolution for D=1. 

- Problem3.py

This code solves the conservative form of hydro equations (1-D Hydro solver) with reflective boundary conditions. It uses the donor cell advection scheme to follow the motion of sound waves in a gas with uniform density and no gravity. The motion starts with a small Gaussian perturbation in density.

To solve the hydro equations we represent the density as "f1" and density times bulk velocity as "f2". Then, the output of this code is one plot showing the evolution of f1 (left panel) where we see a shock developing.

------------------------------------
4. Written explanations
------------------------------------

Problem3.py

1) What happens as you increase the amplitude of the perturbation? Do you see a shock? 

Answer: Increasing the amplitude of the gaussian perturbation results in a right-moving blob that is larger (in amplitude), steeper and wider. The steep portion of the blob looks like the shock front, so yes, I see a shock.

2) If so, what do you think is setting the width of the shock?

Answer: According to my answer above, I see a wider shock when increasing the amplitude of the perturbation. In class we discussed that the shock width is set by the bulk velocity and viscosity. We do not include viscosity explicitly, but there is a numerical viscosity that arises from discretizing the advection terms. It is easier to see though, that the bulk velocity array directly depends on f1 and f2, which are set by the gaussian perturbation. Therefore increasing the amplitude impacts the bulk velocity, which impacts the width of the shock.

------------------------------------
5. Collaborators: Alice Curtin for Problem 3
------------------------------------
