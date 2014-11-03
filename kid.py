#!/usr/bin/python

import cffi
ffi = cffi.FFI()
lib = ffi.dlopen('KiD_1D.so')
ffd.cdef("void save_ptr(char*,void(*)());")
ffi.cdef("void __main_MOD_main_loop();")
lib.__main_MOD_main_loop()
