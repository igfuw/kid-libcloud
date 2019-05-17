import numpy as np
import cffi
from setup import params
import libcloudphxx as libcl
import math
import pdb


ffi = cffi.FFI()
flib = ffi.dlopen('KiD_1D.so')

# Fortran functions (_sp_ means single precision)
ffi.cdef("void __diagnostics_MOD_save_dg_1d_sp_c(float*, int, int, char*, int,     char*, int,      int   );")
#                                                field,  nx,  nz,  name,  namelen, units, unitslen, itime 
ffi.cdef("void __diagnostics_MOD_save_bindata_sp_c(float*, int, char*, int,     char*, int      );")
#                                                  field,  nb,  name,  namelen, units, unitslen                                                                                        
ffi.cdef("void __diagnostics_MOD_save_dg_1d_bin_sp_c(float*, int, int, int, char*, int, char*, int,      int   );")
#                                                    field,  nb,  nx,  nz,  name,  namelen,units, unitslen, itime

ffi.cdef("void __diagnostics_MOD_save_dg_scalar_sp_c(float, char*, int,     char*, int,      int );")
#                                                    scalar, name, namelen, units, unitslen, itime
#    call save_dg(time, 'time', i_dgtime,  units='s',dim='time')

def save_helper(arr):
  # astype() takes keywords arguments for newer numpy versions (1.7?)                       
  try:
    arr = arr.astype(np.float32, copy=False)
  except TypeError:
    arr = arr.astype(np.float32)
  arr_ptr = ffi.cast("float*", arr.__array_interface__['data'][0])
  return arr, arr_ptr
#  # remove the subterrenean level z=0 from output
#  print arr.ndim, arr
#  if (arr.ndim == 2):
#    arr_wo_subterr = np.delete(arr,0,1)
#    print arr_wo_subterr
#    arr_ptr = ffi.cast("float*", arr_wo_subterr.__array_interface__['data'][0])
#    return arr_wo_subterr, arr_ptr
#  else:
#    arr_ptr = ffi.cast("float*", arr.__array_interface__['data'][0])
#    return arr, arr_ptr

def save_dg_scalar(scal, it, name, units):
#  scal = scal.astype(np.float32)
#  scal_cast = ffi.cast("float", scal)
  flib.__diagnostics_MOD_save_dg_scalar_sp_c(
    scal,
    name, len(name),
    units, len(units),
    it
  )

def save_dg(arr, it, name, units):
  arr, arr_ptr = save_helper(arr)
  if (arr.ndim == 2):
    flib.__diagnostics_MOD_save_dg_1d_sp_c(
      arr_ptr, arr.shape[0], arr.shape[1],
      name, len(name),
      units, len(units),
      it
    )
  elif (arr.ndim == 3):
    flib.__diagnostics_MOD_save_dg_1d_bin_sp_c(
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


def diagnostics(particles, arrays, prev_val, it, size_x, size_z, first_timestep):

  # super-droplet concentration per grid cell                               
#  particles.diag_all() 
#  particles.diag_sd_conc()
#  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z), it, "number_of_SDs", "1")

  # recording puddle
  puddle = particles.diag_puddle();
  #save_dg_scalar(puddle[8] * (4./3) * math.pi, it, "accumulated surface precipitation volume", "m^3")
  save_dg_scalar(puddle[8] * (4./3) * math.pi * 1000 / 10, it, "surface precipitation rate", "mm/s") # hardcoded dt_output = 10s

  if first_timestep:
    # temporary arrays (allocating only once)                 
    arrays["tmp_xz"] = np.empty((size_x-2, size_z))
    arrays["mom_0"] = np.empty((params["n_bins"], size_x-2, size_z))
    arrays["mom_3"] = np.empty((params["n_bins"], size_x-2, size_z))

    # upper diameter of a bin with values set using mass-doubling scheme    
    arrays["bins_D_upper"] = (2**np.arange(params["n_bins"]))**(1./3) * params["bin0_D_upper"]
    save_bindata(arrays["bins_D_upper"] * 1e6, "bins_D_upper", "microns")
    save_bindata((arrays["bins_D_upper"]/2)**3 * libcl.common.rho_w * (4./3) * math.pi, "bins_mass_upper", "kg")
    save_bindata(np.diff(np.concatenate([np.zeros(1), arrays["bins_D_upper"]])) * 1e6, "dD", "microns")

    # inferring rain water range as all bigger than cloud water
    params["bins_qr_r20um"] = np.arange(params["bins_qc_r20um"][-1]+1, params["n_bins"])
    #params["bins_qr_r32um"] = np.arange(params["bins_qc_r32um"][-1]+1, params["n_bins"])
    
    
  # T according to the formula used within the library
  #for i in range(0, size_x-2):
  #  for j in range(0, size_z):
  #    arrays["tmp_xz"][i,j] = libcl.common.T(arrays["thetad"][i,j], arrays["rhod"][j])
  #particles.diag_all() 
  #particles.diag_temperature()
  #save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z), it, "T_lib_post_cond", "K")
#  save_dg(arrays["tmp_xz"], it, "T_lib_post_cond", "K")
  #save_dg(arrays["T_lib_ante_cond"], it, "T_lib_ante_cond", "K")

  # RH from the library pre cond
  #save_dg(arrays["RH_lib_ante_cond"], it, "RH_lib_ante_cond", "%")
  # RH from the library pre cond
  #particles.diag_all()
  #particles.diag_RH()
  #save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * 100, it, "RH_lib_post_cond", "%")
  # RH from the library formula with rv post cond and T pre cond (should be the same as RH saved in KiD)
  #for i in range(0, size_x-2):
  #  for j in range(0, size_z):
  #    arrays["tmp_xz"][i,j] = arrays["rhod"][j] * arrays["qv"][i,j] * libcl.common.R_v * arrays["T_lib_ante_cond"][i,j] / libcl.common.p_vs(arrays["T_lib_ante_cond"][i,j]) * 100
  #save_dg(arrays["tmp_xz"], it, "RH_lib_formula_TAnteCond_Rvpostcond", "K")
  #save_dg(arrays["qv"], it, "qv_for_RH_lib_formula_TAnteCond_Rvpostcond_calc", "")

  # saturation vapour pressure from the library formula ante cond
  #save_dg(arrays["psat_lib_formula_ante_cond"], it, "psat_lib_formula_ante_cond", "")

  # pressure from the lib
  #particles.diag_all()
  #particles.diag_pressure()
  #save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * 100, it, "pressure_lib_post_cond", "%")
  #save_dg(arrays["pressure_lib_ante_cond"], it, "pressure_lib_ante_cond", "%")

  # aerosol concentration
  #assert params["bins_qc_r20um"][0] == params["bins_qc_r32um"][0]
  #particles.diag_wet_rng(0, arrays["bins_D_upper"][params["bins_qc_r20um"][0]] / 2)
  #particles.diag_wet_mom(0)
  #save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z), it, "aerosol_number", "/kg") 

  #particles.diag_all()
  #particles.diag_wet_mom(0)
  #liquid_drops_number = np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z)
  #save_dg(liquid_drops_number, it, "liquid_drops_number", "/kg") 
  #particles.diag_all()
  #particles.diag_wet_mom(1)
  #save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * 1e6 / liquid_drops_number, it, "liquid_drops_mean_r", "um") 

  particles.diag_all()
  particles.diag_wet_mom(3)
  liquid_drops_mass = np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  #save_dg(liquid_drops_mass, it, "liquid_drops_mass", "kg/kg")
  save_dg_scalar(np.sum(liquid_drops_mass)*25, it, "liquid water path", "kg/m^2")  # WARNING: hardcoded dz = 25...

  particles.diag_all() 
  particles.diag_accr20()
  accr20 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = accr20.sum()
  save_dg_scalar((new_val - prev_val["accr20"]) / 10, it, "total accretion rate r20", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["accr20"] = new_val

  particles.diag_all() 
  particles.diag_accr25()
  accr25 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = accr25.sum()
  save_dg_scalar((new_val - prev_val["accr25"]) / 10, it, "total accretion rate r25", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["accr25"] = new_val

  particles.diag_all() 
  particles.diag_accr32()
  accr32 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = accr32.sum()
  save_dg_scalar((new_val - prev_val["accr32"]) / 10, it, "total accretion rate r32", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["accr32"] = new_val

  particles.diag_all() 
  particles.diag_acnv20()
  acnv20 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = acnv20.sum()
  save_dg_scalar((new_val - prev_val["acnv20"]) / 10, it, "total autoconversion rate r20", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["acnv20"] = new_val

  particles.diag_all() 
  particles.diag_acnv25()
  acnv25 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = acnv25.sum()
  save_dg_scalar((new_val - prev_val["acnv25"]) / 10, it, "total autoconversion rate r25", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["acnv25"] = new_val

  particles.diag_all() 
  particles.diag_acnv32()
  acnv32 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = acnv32.sum()
  save_dg_scalar((new_val - prev_val["acnv32"]) / 10, it, "total autoconversion rate r32", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["acnv32"] = new_val

  particles.diag_all() 
  particles.diag_revp20()
  revp20 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = revp20.sum()
  save_dg_scalar((new_val - prev_val["revp20"]) / 10, it, "total rain evaporation rate r20", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["revp20"] = new_val

  particles.diag_all() 
  particles.diag_revp25()
  revp25 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = revp25.sum()
  save_dg_scalar((new_val - prev_val["revp25"]) / 10, it, "total rain evaporation rate r25", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["revp25"] = new_val

  particles.diag_all() 
  particles.diag_revp32()
  revp32 =  np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi
  new_val = revp32.sum()
  save_dg_scalar((new_val - prev_val["revp32"]) / 10, it, "total rain evaporation rate r32", "kg/s") #WARNING: hardcoded diag output = 10s
  prev_val["revp32"] = new_val

#  particles.diag_all() 
#  particles.diag_accr32()
#  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "accumulated_accr32_mass", "kg")
#
#  particles.diag_all() 
#  particles.diag_acnv20()
#  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "accumulated_acnv20_mass", "kg")
#
#  particles.diag_all() 
#  particles.diag_acnv32()
#  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "accumulated_acnv32_mass", "kg")

  # binned wet spectrum                                      
  r_min = 0.
  for i in range(params["n_bins"]):
    # selecting range                                 
    r_max = arrays["bins_D_upper"][i] / 2
    particles.diag_wet_rng(r_min, r_max)
    r_min = r_max
    # computing 1-st moment                                  
    particles.diag_wet_mom(0)
    arrays["mom_0"][i,:,:] = np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z)
    # computing 3-rd moment                          
    particles.diag_wet_mom(3)
    arrays["mom_3"][i,:,:] = np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z)

  # saving binned concentrations
  save_dg(arrays["mom_0"], it, "cloud_bin_number", "/kg")

  # saving summed concentrations
  #save_dg(np.sum(arrays["mom_0"][params["bins_qc_r20um"],:,:], axis=0), it, "cloud_number_r20um", "/kg")
  #save_dg(np.sum(arrays["mom_0"][params["bins_qc_r32um"],:,:], axis=0), it, "cloud_number_r32um", "/kg")
  #save_dg(np.sum(arrays["mom_0"][params["bins_qr_r20um"],:,:], axis=0), it, "rain_number_r20um", "/kg")
  #save_dg(np.sum(arrays["mom_0"][params["bins_qr_r32um"],:,:], axis=0), it, "rain_number_r32um", "/kg")

  # saving binned masses
  arrays["mom_3"] *= libcl.common.rho_w * (4./3) * math.pi
  save_dg(arrays["mom_3"], it, "cloud_bin_mass", "kg/kg")
  
  # saving summed masses
  #save_dg(np.sum(arrays["mom_3"][params["bins_qc_r20um"],:,:], axis=0), it, "cloud_mass_r20um", "kg/kg")
  #save_dg(np.sum(arrays["mom_3"][params["bins_qc_r32um"],:,:], axis=0), it, "cloud_mass_r32um", "kg/kg")
  #save_dg(np.sum(arrays["mom_3"][params["bins_qr_r20um"],:,:], axis=0), it, "rain_mass_r20um", "kg/kg")
  #save_dg(np.sum(arrays["mom_3"][params["bins_qr_r32um"],:,:], axis=0), it, "rain_mass_r32um", "kg/kg")

  #save_dg_scalar(np.sum(arrays["mom_3"][params["bins_qr_r20um"],:,:]) * 25, it, "rain water path r20 bins", "kg/kg")
  #save_dg_scalar(np.sum(arrays["mom_3"][params["bins_qr_r25um"],:,:]) * 25, it, "rain water path r25", "kg/kg")
  #save_dg_scalar(np.sum(arrays["mom_3"][params["bins_qr_r32um"],:,:]) * 25, it, "rain water path r32", "kg/kg")

  particles.diag_wet_rng(20e-6, 1)
  particles.diag_wet_mom(3)
  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z), it, "rain_water_r20", "kg/kg") 
  save_dg_scalar(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z).sum() * libcl.common.rho_w * (4./3) * math.pi * 25, it, "rain_water_path_r20", "kg/m^2") 
  particles.diag_wet_rng(25e-6, 1)
  particles.diag_wet_mom(3)
  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "rain_water_r25", "kg/kg") 
  save_dg_scalar(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z).sum() * libcl.common.rho_w * (4./3) * math.pi * 25, it, "rain_water_path_r25", "kg/m^2") 
  particles.diag_wet_rng(32e-6, 1)
  particles.diag_wet_mom(3)
  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "rain_water_r32", "kg/kg") 
  save_dg_scalar(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z).sum() * libcl.common.rho_w * (4./3) * math.pi * 25, it, "rain_water_path_r32", "kg/m^2") 


  particles.diag_wet_rng(0.5e-6, 20e-6)
  particles.diag_wet_mom(3)
  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "cloud_water_r20", "kg/kg") 
  particles.diag_wet_rng(0.5e-6, 25e-6)
  particles.diag_wet_mom(3)
  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "cloud_water_r25", "kg/kg") 
  particles.diag_wet_rng(0.5e-6, 32e-6)
  particles.diag_wet_mom(3)
  save_dg(np.frombuffer(particles.outbuf()).copy().reshape(size_x-2, size_z) * libcl.common.rho_w * (4./3) * math.pi, it, "cloud_water_r32", "kg/kg") 

  # binned dry spectrum? - TODO                       
  # ...                                                               
