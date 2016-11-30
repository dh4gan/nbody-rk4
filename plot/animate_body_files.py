# Written 30/11/16 by dh4gan
# Reads in data from nbody_rk4 (individual body files)
# Animates the result

import io_nbody_rk4 as io


prefix = raw_input('What is the run prefix?')
prefix = prefix+'*'
io.animate_bodies(prefix)