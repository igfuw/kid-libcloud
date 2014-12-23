import numpy as np
import cffi
from setup import params
import libcloudphxx as libcl
import math
import pdb

ffi = cffi.FFI()
flib = ffi.dlopen('KiD_SC_2D.so')

# Fortran functions (_sp_ means single precision)
ffi.cdef("void __diagnostics_MOD_save_dg_2d_sp_c(float*, int, int, char*, int,     char*, int,      int   );")
#                                                field,  nx,  nz,  name,  namelen, units, unitslen, itime 
ffi.cdef("void __diagnostics_MOD_save_bindata_sp_c(float*, int, char*, int,     char*, int      );")
#                                                  field,  nb,  name,  namelen, units, unitslen                                                                                        
ffi.cdef("void __diagnostics_MOD_save_dg_2d_bin_sp_c(float*, int, int, int, char*, int, char*, int,      int   );")
#                                                    field,  nb,  nx,  nz,  name,  namelen,units, unitslen, itime

def save_helper(arr):
  # astype() takes keywords arguments for newer numpy versions (1.7?)                       
  try:
    arr = arr.astype(np.float32, copy=False)
  except TypeError:
    arr = arr.astype(np.float32)
  arr_ptr = ffi.cast("float*", arr.__array_interface__['data'][0])
  return arr, arr_ptr


def save_dg(arr, it, name, units):
  arr, arr_ptr = save_helper(arr)
  if (arr.ndim == 2):
    flib.__diagnostics_MOD_save_dg_2d_sp_c(
      arr_ptr, arr.shape[0], arr.shape[1],
      name, len(name),
      units, len(units),
      it
    )
  elif (arr.ndim == 3):
    flib.__diagnostics_MOD_save_dg_2d_bin_sp_c(
      arr_ptr, arr.shape[0], arr.shape[1], arr.shape[2],
      name, len(name),
      units, len(units),
      it
    )
  else:
    assert(False)


def save_bindata(arr, name, unit):
  assert(arr.ndim == 1)
  arr, arr_ptr = save_helper(arr)
  flib.__diagnostics_MOD_save_bindata_sp_c(
    arr_ptr, arr.shape[0],
    name, len(name),
    unit, len(unit)
  )


def diagnostics(particles, arrays, it, size_x, size_z, first_timestep):

  tmp = np.empty((size_x-2, size_z))
  tmp[:,:] = arrays["qv"]
  save_dg(tmp, it, "aqq", "J")

  # super-droplet concentration per grid cell                               
  # TODO: select all particles?                                 
  particles.diag_sd_conc()
  save_dg(np.frombuffer(particles.outbuf()).reshape(size_x-2, size_z), it, "sd_conc", "1")

  # temporary arrays (allocating only once)                 
  if first_timestep:
    arrays["tmp_xz"] = np.empty((size_x-2, size_z))
    arrays["mom_0"] = np.empty((params["n_bins"], size_x-2, size_z))
    arrays["mom_3"] = np.empty((params["n_bins"], size_x-2, size_z))
    # upper diameter of a bin with values set using mass-doubling scheme    
    arrays["bins_D_upper"] = (2**np.arange(params["n_bins"]))**(1./3) * params["bin0_D_upper"]
    save_bindata(arrays["bins_D_upper"] * 1e6, "bins_D_upper", "microns")
    save_bindata((arrays["bins_D_upper"]/2)**3 * libcl.common.rho_w * (4./3) * math.pi, "bins_mass_upper", "kg")
    save_bindata(np.diff(np.concatenate([np.zeros(1), arrays["bins_D_upper"]])) * 1e6, "dD", "microns")
    
  # T according to the formula used within the library
  for i in range(0, size_x-2):
    for j in range(0, size_z):
      arrays["tmp_xz"][i,j] = libcl.common.T(arrays["thetad"][i,j], arrays["rhod"][j])
  save_dg(arrays["tmp_xz"], it, "T_lib_post_cond", "K")

  # RH according to the formula used within the library
  for i in range(0, size_x-2):
    for j in range(0, size_z):
      arrays["tmp_xz"][i,j] = arrays["rhod"][j] * arrays["qv"][i,j] * libcl.common.R_v * arrays["tmp_xz"][i,j] / libcl.common.p_vs(arrays["tmp_xz"][i,j])
  save_dg(arrays["tmp_xz"], it, "RH_lib_post_cond", "K")

  # binned wet spectrum                                      
  r_min = 0.
  for i in range(params["n_bins"]):
    # selecting range                                 
    r_max = arrays["bins_D_upper"][i] / 2
    particles.diag_wet_rng(r_min, r_max)
    r_min = r_max
    # computing 1-st moment                                  
    particles.diag_wet_mom(0)
    arrays["mom_0"][i,:,:] = np.frombuffer(particles.outbuf()).reshape(size_x-2, size_z)
    # computing 3-rd moment                          
    particles.diag_wet_mom(3)
    arrays["mom_3"][i,:,:] = np.frombuffer(particles.outbuf()).reshape(size_x-2, size_z)

  save_dg(arrays["mom_0"], it, "cloud_bin_number", "/kg")
  arrays["mom_3"] *= libcl.common.rho_w * (4./3) * math.pi
  save_dg(arrays["mom_3"], it, "cloud_bin_mass", "kg/kg")

  # binned dry spectrum? - TODO                       
  # ...                                                               
