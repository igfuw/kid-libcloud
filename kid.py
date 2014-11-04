#!/usr/bin/python

import numpy as np
import cffi
ffi = cffi.FFI()

real_t = np.float64

@ffi.callback("void(double*, int, int)")
def hello(tablica, nz, nx):
  print "hello from Python", tablica, nx, nz
  array = np.frombuffer(ffi.buffer(tablica, nx*nz*np.dtype(real_t).itemsize), dtype=real_t)
  print array.reshape((nx,nz))[1:-1,:]

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# Fortran functions
ffi.cdef("void __main_MOD_main_loop();")

lib = ffi.dlopen('KiD_1D.so')

# storing pointers to Python functions
lib.save_ptr("hello.ptr", hello)

# running Fortran stuff
lib.__main_MOD_main_loop()
