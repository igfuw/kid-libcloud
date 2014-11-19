#!/usr/bin/python

import numpy as np
import cffi
ffi = cffi.FFI()

real_t = np.float64

@ffi.callback("void(int, int, int, double*, double*, double*, double*, double*)")
def micro_step(it_diag, size_z, size_x, th_ar, qv_ar, rho_ar, uh_ar, wh_ar):
  print "micro_step from Python", th_ar, size_z, size_x
  print "i_dgtime z fortrana", it_diag
  array = np.frombuffer(ffi.buffer(th_ar, size_x*size_z*np.dtype(real_t).itemsize), dtype=real_t)
  #print array.reshape((size_x,size_z))[1:-1,:]

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# Fortran functions
ffi.cdef("void __main_MOD_main_loop();")

lib = ffi.dlopen('KiD_SC_2D.so')

# storing pointers to Python functions
lib.save_ptr("/tmp/micro_step.ptr", micro_step)

# running Fortran stuff
# note: not using command line arguments, namelist name hardcoded in
#       kid_a_setup/namelists/SC_2D_input.nml 
lib.__main_MOD_main_loop()
