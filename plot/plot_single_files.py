# Written 29/11/16 by dh4gan
# Plots data from nbody_rk4
# Handles single body file formats


import io_nbody_rk4 as io
import filefinder as ff


prefix = raw_input('What is the run prefix? ')

filename = ff.find_local_input_files(prefix+'*')
io.plot_body(filename, 'single')
