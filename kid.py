#!/usr/bin/python

import numpy as np
import cffi
import libcloudphxx as libcl
from libcloudphxx.common import R_v, R_d, c_pd
import pdb

# CFFI stuff
ffi = cffi.FFI()
lib = ffi.dlopen('KiD_SC_2D.so')

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# Fortran functions
ffi.cdef("void __main_MOD_main_loop();")
ffi.cdef("void __diagnostics_MOD_save_dg_2d_sp_c(float*, int, int, char*, int, char*, int, int);")

# object storing super-droplet model state (to be initialised)
prtcls = False

# dictionary of simulation parameters
params = {
  "real_t" : np.float64,
  "backend" : libcl.lgrngn.backend_t.serial,
  "sd_conc" : 128.,
  "kappa" : .61,
  "meanr" : .04e-6,
  "gstdv" : 1.4,
  "n_tot" : 100e6
}

arrays = {}
dt, dx, dz = 0, 0, 0
first_timestep = True
last_diag = -1
opts = libcl.lgrngn.opts_t()
opts.sstp_cond = 1
opts.sstp_coal = 1
opts.cond = True
opts.coal = True
opts.adve = True
opts.sedi = True
opts.chem = False
opts.kernel = libcl.lgrngn.kernel_t.geometric

def lognormal(lnr):
  from math import exp, log, sqrt, pi
  return params["n_tot"] * exp(
    -pow((lnr - log(params["meanr"])), 2) / 2 / pow(log(params["gstdv"]),2)
  ) / log(params["gstdv"]) / sqrt(2*pi);

def ptr2np(ptr, size_x, size_z):
  numpy_ar = np.frombuffer(
    ffi.buffer(ptr, size_x*size_z*np.dtype(params["real_t"]).itemsize),
    dtype=params["real_t"]
  ).reshape(size_x, size_z)
  return numpy_ar.squeeze()

def th_kid2dry(th, rv):
  return th * (1 + rv * R_v / R_d)**(R_d/c_pd)

def th_dry2kid(th_d, rv):
  return th_d * (1 + rv * R_v / R_d)**(-R_d/c_pd)

def rho_kid2dry(rho, rv):
  return rho / (1 + rv) #TODO: I'm assuming that KiD uses rho

def save_dg(arr, it, name, units):
  # astype() takes keywords arguments for newer numpy versions (1.7?)
  try: 
    arr = arr.astype(np.float32, copy=False) 
  except TypeError:
    arr = arr.astype(np.float32)
  arr_ptr = ffi.cast("float*", arr.__array_interface__['data'][0])
  lib.__diagnostics_MOD_save_dg_2d_sp_c(
    arr_ptr, arr.shape[0], arr.shape[1], 
    name, len(name), 
    units, len(units),
    it
  )

def diagnostics(particles, it, size_x, size_z):
  tmp = np.empty((size_x-2, size_z))
  tmp[:,:] = arrays["qv"]
  save_dg(tmp, it, "aqq", "J")

  # TODO: select all particles
  particles.diag_sd_conc()
  save_dg(np.frombuffer(particles.outbuf()).reshape(size_x-2, size_z), it, "sd_conc", "1")
  

@ffi.callback("void(int, float, int, int, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*)")
def micro_step(it_diag, dt, size_z, size_x, th_ar, qv_ar, rhof_ar, rhoh_ar, 
               uf_ar, uh_ar, wf_ar, wh_ar, xf_ar, zf_ar, xh_ar, zh_ar, 
               tend_th_ar, tend_qv_ar, err):
  # global should be used for all variables defined in "if first_timestep"  
  global prtcls, dx, dz, first_timestep, last_diag

  print " w pythonie na poczatku err", err
  err_py = ffi.cast("int *",err) 
  try:
    # superdroplets: initialisation (done only once)
    if first_timestep:

      arrx = ptr2np(xf_ar, size_x, 1)
      arrz = ptr2np(zf_ar, 1, size_z)

      # checking if grids are equal
      np.testing.assert_almost_equal((arrx[1:]-arrx[:-1]).max(), (arrx[1:]-arrx[:-1]).min(), decimal=7)
      np.testing.assert_almost_equal((arrz[1:]-arrz[:-1]).max(), (arrz[1:]-arrz[:-1]).min(), decimal=7)
      dx = arrx[1] - arrx[0]                            
      dz = arrz[1] - arrz[0]                            
      
      opts_init = libcl.lgrngn.opts_init_t()
      opts_init.dt = dt
      opts_init.nx, opts_init.nz = size_x - 2, size_z
      opts_init.dx, opts_init.dz = dx, dz 
      opts_init.x1, opts_init.z1 = dx * opts_init.nx, dz * opts_init.nz
      opts_init.sd_conc_mean = params["sd_conc"]
      opts_init.dry_distros = { params["kappa"] : lognormal }
      
      prtcls = libcl.lgrngn.factory(params["backend"], opts_init)
      
      # allocating arrays for those variables that are not ready to use
      # (i.e. either different size or value conversion needed)
      for name in ("thetad", "qv"):
        arrays[name] = np.empty((size_x-2, size_z))
      arrays["rhod"] = np.empty((size_x-2, size_z)) #TODO: add to loop or change back to one dim
      arrays["rhod_Cx"] = np.empty((size_x-1, size_z))
      arrays["rhod_Cz"] = np.empty((size_x-2, size_z+1))
#end if
    
    # mapping local NumPy arrays to the Fortran data locations   
    arrays["qv"][:,:] = ptr2np(qv_ar, size_x, size_z)[1:-1, :]
    arrays["thetad"][:,:] = th_kid2dry(ptr2np(th_ar, size_x, size_z)[1:-1, :], arrays["qv"][:,:])
    arrays["rhod"][:,:] = rho_kid2dry(ptr2np(rhof_ar, 1, size_z)[:], arrays["qv"][:,:])
    
    arrays["rhod_Cx"][:,:] = ptr2np(uh_ar, size_x, size_z)[:-1, :]
    assert (arrays["rhod_Cx"][0,:] == arrays["rhod_Cx"][-1,:]).all()
    arrays["rhod_Cx"] *= arrays["rhod"][0] * dt / dx #TODO - rhod should be in uh_ar places??
    
    arrays["rhod_Cz"][:, 1:] = ptr2np(wh_ar, size_x, size_z)[1:-1, :] 
    arrays["rhod_Cz"][:, 0 ] = 0
    arrays["rhod_Cz"][:, 1:] *= rho_kid2dry(ptr2np(rhoh_ar, 1, size_z), arrays["qv"][:,:]) * dt / dz

  
    if first_timestep:
      prtcls.init(arrays["thetad"], arrays["qv"], arrays["rhod"], arrays["rhod_Cx"], arrays["rhod_Cz"]) 
      diagnostics(prtcls, 1, size_x, size_z) # writing down state at t=0

    # superdroplets: all what have to be done within a timestep
    prtcls.step_sync(opts, arrays["thetad"], arrays["qv"],  arrays["rhod"]) #TODO: courants...
    prtcls.step_async(opts)

    # calculating tendency for theta (first converting back to non-dry theta
    ptr2np(tend_th_ar, size_x, size_z)[1:-1, :] = (
      ptr2np(th_ar, size_x, size_z)[1:-1, :] -   # old
      th_dry2kid(arrays["thetad"], arrays["qv"]) # new
      ) / dt #TODO: check if dt needed
    
    # calculating tendency for qv
    ptr2np(tend_qv_ar, size_x, size_z)[1:-1, :] = (
      ptr2np(qv_ar, size_x, size_z)[1:-1, :] - # old                
      arrays["qv"]                             # new 
      ) / dt #TODO: check if dt needed    

    # diagnostics
    if last_diag < it_diag:
      diagnostics(prtcls, it_diag, size_x, size_z)
      last_diag = it_diag

    first_timestep = False
    err_py[0] = 0
  except:
    print " w except"
    err_py[0] = 1
  
  
# storing pointers to Python functions
lib.save_ptr("/tmp/micro_step.ptr", micro_step)

# running Fortran stuff
# note: not using command line arguments, namelist name hardcoded in
#       kid_a_setup/namelists/SC_2D_input.nml 
lib.__main_MOD_main_loop()
