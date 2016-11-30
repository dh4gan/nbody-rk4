#####################################################
###                                               ###
###           Makefile for nbody_rk4              ###
###                                               ###
###         Duncan H. Forgan 18/11/2016           ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables (standard):
FC     = gfortran 

# For files generated from stacpolly use these flags
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -Wall -fbounds-check

# Create object files:
%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = nbodymodule.f90 main.f90 calc_system_properties.f90 endrun.f90 \
	initial.f90 integrate.f90 calc_acceleration.f90 calc_drag_terms.f90 \
	calc_grav_acceleration.f90 \
	orbits.f90 output.f90 timestep.f90
OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: nbody_rk4

nbody_rk4:  $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean: 
	\rm *.o *.mod nbody_rk4

# End Makefile
