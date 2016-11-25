### N-Body integration via 4th-order Runge Kutta ###
----------------------------------------------------

This FORTRAN 90 repository sets up N bodies and solves their equation
of motion numerically via the 4th-order Runge Kutta method (RK4) with a global (adaptive) timestep.  

Note - this code was developed so that the algorithms contained within could be integrated into another code.

It was written standalone as a means of testing the above algorithms.

I've done my best to make this simple to use and user-facing, but I don't pretend that it's entirely user-friendly.

The code also calculates orbital elements (by default around the centre of mass).  The principal force is
Newtonian Gravitation, but drag terms are also implementable.  These terms mimic the drag forces experienced
by particles in gaseous discs (radial inward migration, eccentricity and inclination damping).

## Outputs ##
The code can produce two types of output:

1) One file per particle, to easily plot particle properties vs time
2) Snapshot files, where all particle data at a single timestep is recorded


## Compilation and Execution ##

This code was developed and tested on gfortran, and compiled with standard Makefile software.

The code is compiled using the `Makefile` contained within the repo, and run with the command

`> ./nbody_rk4 <paramfile> `

Where <paramfile> is the parameters file needed to setup the N Body system (example: `nbody_rk4.params`)



