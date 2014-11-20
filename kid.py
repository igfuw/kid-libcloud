#!/usr/bin/python

import numpy as np
import cffi
import libcloudphxx as libcl
import pdb

ffi = cffi.FFI()

# object storing super-droplet model state (to be initialised)
prtcls = False

# dictionary of simulation parameters
params = {
  "real_t" : np.float64,
  "backend" : libcl.lgrngn.backend_t.serial,
  "sd_conc" : 128,
  "kappa" : .61,
  "meanr" : .04e-6,
  "gstdv" : 1.4,
  "n_tot" : 100e6
}

def lognormal(lnr):
  from math import exp, log, sqrt, pi
  print lnr
  return params["n_tot"] * exp(
    -pow((lnr - log(params["meanr"])), 2) / 2 / pow(log(params["gstdv"]),2)
  ) / log(params["gstdv"]) / sqrt(2*pi);

def ptr2np(ptr, size_x, size_z):
  # TODO: halo!, strides                                                                
  numpy_ar = np.frombuffer(
    ffi.buffer(ptr, size_x*size_z*np.dtype(params["real_t"]).itemsize),
    dtype=params["real_t"]
  ).reshape(size_x, size_z)
  return numpy_ar



@ffi.callback("void(int, int, int, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*)")
def micro_step(it_diag, size_z, size_x, th_ar, qv_ar, rho_ar, 
               uf_ar, uh_ar, wf_ar, wh_ar, xf_ar, zf_ar, xh_ar, zh_ar):
  global prtcls

  
  print "in python::micro_step() from Python"#, prtcls #, th_ar, size_z, size_x

  # superdroplets: initialisation (done only once)
  if not prtcls:
    print "initialisation!"
    opts_init = libcl.lgrngn.opts_init_t()

    opts_init.sd_conc = params["sd_conc"]
    opts_init.dry_distros = { params["kappa"] : lognormal }

    prtcls = libcl.lgrngn.factory(params["backend"], opts_init)
    #prtcls.init()

  # mapping local NumPy arrays to the Fortran data locations   
  #print " x, x_hlf", ptr2np(xf_ar, size_x, 1), ptr2np(xh_ar, size_x, 1)
  #print "z, z_half", ptr2np(zf_ar, 1, size_z), ptr2np(zh_ar, 1, size_z)
  #print "w, w_half", ptr2np(wf_ar, size_x, size_z), ptr2np(wh_ar, size_x, size_z)


  # superdroplets: all what have to be done within a timestep
  #prtcls.step_sync()
  #prtcls.step_async()

  #print "i_dgtime z fortrana", it_diag
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
