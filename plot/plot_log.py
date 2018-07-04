# Written 29/11/16 by dh4gan
# Plots data from nbody_rk4
# Handles log body file formats


import io_nbody_rk4 as io
import filefinder as ff

filename = ff.find_local_input_files('*.log')
io.plot_body(filename, 'log')
