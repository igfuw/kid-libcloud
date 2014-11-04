#!/usr/bin/python

import cffi
ffi = cffi.FFI()

@ffi.callback("void(double*)")
def hello(tablica):
  print "hello from Python", tablica

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# Fortran functions
ffi.cdef("void __main_MOD_main_loop();")

lib = ffi.dlopen('KiD_1D.so')

# storing pointers to Python functions
lib.save_ptr("hello.ptr", hello)

# running Fortran stuff
lib.__main_MOD_main_loop()
