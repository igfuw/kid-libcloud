#!/usr/bin/python

import numpy as np
import cffi
import traceback
import libcloudphxx as libcl
from libcloudphxx.common import R_v, R_d, c_pd, eps, p_1000
from setup import params, opts
import diagnostics_1D as dg
import os
import json
import pdb

from argparse import ArgumentParser


ptrfname = "/tmp/micro_step-" + str(os.getuid()) + "-" + str(os.getpid()) + ".ptr"

# CFFI stuff
ffi = cffi.FFI()
flib = ffi.dlopen('KiD_1D.so')
clib = ffi.dlopen('ptrutil.so')

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# Fortran functions (_sp_ means single precision)
ffi.cdef("void __main_MOD_main_loop();")

# object storing super-droplet model state (to be initialised)
prtcls = False

arrays = {}
timestep = 0
last_diag = -1

#parser for overriding values from setup.py with command-line arguments
prsr = ArgumentParser(add_help=True, description='1D kidA case')
prsr.add_argument('--n_tot', required=False, type=float, default=params["n_tot"], help='initial aerosol concentation at STP [1/kg_dry_air]')
prsr.add_argument('--spinup_rain', required=False, type=float, default=params["spinup_rain"], help='time, after which coalescence and sedimentation are turned on [s]')
args = prsr.parse_args()

params["n_tot"] = args.n_tot # * 1.225 / 1. # 1.225 is air density at stp, 1 is the actual density
params["spinup_rain"] = args.spinup_rain
#savings some parameters from setup.py file and libcl revision number
params_write = params.copy()

# converting numpy objects to lists or strings, so json can save them
for key_ar in ["bins_qc_r20um", "bins_qc_r32um"]:
  params_write[key_ar] = params[key_ar].tolist()
for key_str in ["real_t"]:
  params_write[key_str] = str(params[key_str])
params_write["libcloudph_git_rev"] = libcl.git_revision

file_out = open("output/python_setup.txt", "w")
json.dump(params_write, file_out)
file_out.close()

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
  # KiD seems to define rho as (p_v + p_d) / (R_d T)
  return rho / (1 + rv / eps) 

@ffi.callback("bool(int, float, int, int, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*)")
def micro_step(it_diag, dt, size_z, size_x, th_ar, qv_ar, rhof_ar, rhoh_ar, exner_ar, 
               uf_ar, uh_ar, wf_ar, wh_ar, xf_ar, zf_ar, xh_ar, zh_ar, tend_th_ar, tend_qv_ar, rh_ar):
  try:
    # global should be used for all variables defined in "if first_timestep"  
    global prtcls, dx, dz, timestep, last_diag
    #pdb.set_trace()
    # superdroplets: initialisation (done only once)
    if timestep == 0:
      # first, removing the no-longer-needed pointer file
      os.unlink(ptrfname)

      arrx = ptr2np(xf_ar, size_x, 1)
      arrz = ptr2np(zf_ar, 1, size_z)
 
      # checking if grids are equal
      np.testing.assert_almost_equal((arrx[1:]-arrx[:-1]).max(), (arrx[1:]-arrx[:-1]).min(), decimal=7)
      np.testing.assert_almost_equal((arrz[1:]-arrz[:-1]).max(), (arrz[1:]-arrz[:-1]).min(), decimal=7)
      dx = 1.
      dz = arrz[1] - arrz[0]                            
      
      opts_init = libcl.lgrngn.opts_init_t()
      opts_init.dt = dt
      opts_init.nx, opts_init.nz = size_x - 2, size_z
      opts_init.dx, opts_init.dz = dx, dz 
      opts_init.z0 = dz # skipping the first sub-terrain level
      opts_init.x1, opts_init.z1 = dx * opts_init.nx, dz * opts_init.nz
      opts_init.sd_conc = int(params["sd_conc"])
      opts_init.dry_distros = { params["kappa"] : lognormal }
      opts_init.sstp_cond, opts_init.sstp_coal = params["sstp_cond"], params["sstp_coal"]
      opts_init.terminal_velocity = libcl.lgrngn.vt_t.beard77fast
      opts_init.kernel = libcl.lgrngn.kernel_t.hall_davis_no_waals
      opts_init.adve_scheme = libcl.lgrngn.as_t.pred_corr
      opts_init.exact_sstp_cond = 1
      opts_init.sd_conc_large_tail = 1
      opts_init.n_sd_max = int(1.2 * opts_init.nx*opts_init.nz*opts_init.sd_conc)
      opts_init.periodic_z = 0
      opts_init.aerosol_independent_of_rhod = 1 # set to true, because rhod is supposed to be =1, but we cannot pass rhod=1 as it is not in agreement with the values of p and theta and would lead to wrong T,RH,etc...

      print "nx = ", opts_init.nx
      print "nz = ", opts_init.nz
      print "dx = ", opts_init.dx
      print "dz = ", opts_init.dz
      print "dt = ", opts_init.dt
      
      try:
        print("Trying with multi_CUDA backend..."),
        prtcls = libcl.lgrngn.factory(libcl.lgrngn.backend_t.multi_CUDA, opts_init)
        print (" OK!")
      except:
        print (" KO!")
        try:
          print("Trying with CUDA backend..."),
          prtcls = libcl.lgrngn.factory(libcl.lgrngn.backend_t.CUDA, opts_init)
          print (" OK!")
        except:
          print (" KO!")
          try:
            print("Trying with OpenMP backend..."),
            prtcls = libcl.lgrngn.factory(libcl.lgrngn.backend_t.OpenMP, opts_init)
            print (" OK!")
          except:
            print (" KO!")
            print("Trying with serial backend..."),
            prtcls = libcl.lgrngn.factory(libcl.lgrngn.backend_t.serial, opts_init)
            print (" OK!")
    
      # allocating arrays for those variables that are not ready to use
      # (i.e. either different size or value conversion needed)
      for name in ("thetad", "qv", "p_d", "T_kid", "rhod_kid"):
	arrays[name] = np.empty((opts_init.nx, opts_init.nz))
      arrays["rhod"] = np.empty((opts_init.nz,))
      arrays["Cx"] = np.empty((opts_init.nx+1, opts_init.nz))
      arrays["Cz"] = np.empty((opts_init.nx, opts_init.nz+1))
      arrays["RH_lib_ante_cond"] = np.empty((opts_init.nx, opts_init.nz))
      arrays["pressure_lib_ante_cond"] = np.empty((opts_init.nx, opts_init.nz))
      arrays["T_lib_ante_cond"] = np.empty((opts_init.nx, opts_init.nz))
      arrays["psat_lib_formula_ante_cond"] = np.empty((opts_init.nx, opts_init.nz))
      arrays["RH_kid"] = np.empty((opts_init.nx, opts_init.nz))

    # defining qv and thetad (in every timestep) 
    arrays["qv"][:,:] = ptr2np(qv_ar, size_x, size_z)[1:-1, :]
    #arrays["thetad"][:,:] = th_kid2dry(ptr2np(th_ar, size_x, size_z)[1:-1, :], arrays["qv"][:,:])

    # it seems (c.f. src/test_cases.f90:863 which suggest that exner=(p_d / p_0)^(R_d/c_p)) that theta, p and rhod calculated in KiD are actually dry air values
    # so there's no need to make them dry one more time; BTW: this makes the RH calculated in KiD not correct, because
    # it should be calculated using total pressure and not dry pressure
    arrays["thetad"][:,:] = ptr2np(th_ar, size_x, size_z)[1:-1, :]
    arrays["p_d"][:,:] = (ptr2np(exner_ar, size_x, size_z)[1:-1, :])**(c_pd/R_d) * p_1000# * eps / (eps + arrays["qv"])
    arrays["T_kid"][:,:] = ptr2np(exner_ar, size_x, size_z)[1:-1, :] * ptr2np(th_ar, size_x, size_z)[1:-1, :]
    arrays["rhod_kid"] = arrays["p_d"] / arrays["T_kid"] / R_d
    arrays["RH_kid"][:,:] = ptr2np(rh_ar, size_x, size_z)[1:-1, :]
    #pdb.set_trace()

    # finalising initialisation
    # density in KiD (rhof) is overriden to 1, but it is inconsistent with exner pressure and temperature
    # and this leads to wrong RH calculations in libcloud
    # therefore we should use the rhod_kid density and multiply diags by rhod_kid ?!
    if timestep == 0:
      arrays["rhod"][:] = arrays["rhod_kid"][0,:] # set rhod to be in agreement with other profiles and exner function
      #arrays["rhod"][:] = rho_kid2dry(ptr2np(rhof_ar, 1, size_z)[:], arrays["qv"][0,:])  # pass rhod from KiD, i.e. rhod=1

     
    #arrays["Cx"][:,:] = ptr2np(uh_ar, size_x, size_z)[:-1, :] * dt / dx 
    assert (arrays["Cx"][0,:] == arrays["Cx"][-1,:]).all()

    # putting meaningful values to the sub-terain level (to avoid segfault from library)
    arrays["Cz"][:, 0] = 0.
    arrays["qv"][:, 0] = arrays["qv"][:, 1]
    arrays["thetad"][:,0] = arrays["thetad"][:,1]
    arrays["rhod"][0] = arrays["rhod"][1]

    arrays["Cz"][:, 1:] = ptr2np(wh_ar, size_x, size_z)[1:-1, :] * dt / dz

    if timestep == 0:
      prtcls.init(arrays["thetad"], arrays["qv"], arrays["rhod"], Cx = arrays["Cx"], Cz = arrays["Cz"]) 
      dg.diagnostics(prtcls, arrays, 1, size_x, size_z, timestep == 0) # writing down state at t=0

    # spinup period logic
    opts.sedi = opts.coal = timestep >= params["spinup_rain"]
 #   if timestep >= params["spinup_smax"]: opts.RH_max = 44

    # saving RH for the output file
    prtcls.diag_all()
    prtcls.diag_RH()
    arrays["RH_lib_ante_cond"] = np.frombuffer(prtcls.outbuf()).reshape(size_x - 2, size_z) * 100

    prtcls.diag_all()
    prtcls.diag_pressure()
    arrays["pressure_lib_ante_cond"] = np.frombuffer(prtcls.outbuf()).reshape(size_x - 2, size_z) * 100
    
#    if timestep == 0:
#      arrays["rhod"][:] = ptr2np(rhof_ar, 1, size_z)[:]
#      print "rhof_ar: ", arrays["rhod"]
#      arrays["rhod"][:] = rho_kid2dry(ptr2np(rhof_ar, 1, size_z)[:], arrays["qv"][0,:])
#      print "rhod: ", arrays["rhod"]
#      print "rhod_kid: ", arrays["rhod_kid"]



#    for i in range(0, prtcls.opts_init.nx):
#      for j in range(0, prtcls.opts_init.nz):
#        arrays["T_lib_ante_cond"][i,j] = libcl.common.T(arrays["thetad"][i,j], arrays["rhod_kid"][i,j])
#    print "T_lib_ante_cond with rhod_kid: ", arrays["T_lib_ante_cond"]
    for i in range(0, prtcls.opts_init.nx):
      for j in range(0, prtcls.opts_init.nz):
        arrays["T_lib_ante_cond"][i,j] = libcl.common.T(arrays["thetad"][i,j], arrays["rhod"][j])


    for i in range(0, prtcls.opts_init.nx):
      for j in range(0, prtcls.opts_init.nz):
        arrays["psat_lib_formula_ante_cond"][i,j] =  libcl.common.p_vs(arrays["T_lib_ante_cond"][i,j])

    if timestep == 0:
      print "rhod: ", arrays["rhod"]
      print "thetad: ", arrays["thetad"]
      print "qv: ", arrays["qv"]
      print "RH_lib_ante_cond: ", arrays["RH_lib_ante_cond"]
      print "T_lib_ante_cond: ", arrays["T_lib_ante_cond"]
      print "T_kid: ", arrays["T_kid"]

    # superdroplets: all what have to be done within a timestep
    prtcls.step_sync(opts, arrays["thetad"], arrays["qv"], Cx = arrays["Cx"], Cz = arrays["Cz"]) 
    #prtcls.step_sync(opts, arrays["thetad"], arrays["qv"], Cx = arrays["Cx"], Cz = arrays["Cz"], RH = arrays["RH_kid"]) 
    #prtcls.step_sync(opts, arrays["thetad"], arrays["qv"], Cx = arrays["Cx"], Cz = arrays["Cz"], T = arrays["T_kid"]) 
    #prtcls.step_sync(opts, arrays["thetad"], arrays["qv"], Cx = arrays["Cx"], Cz = arrays["Cz"], RH = arrays["RH_kid"], T = arrays["T_kid"]) 

    prtcls.step_async(opts)


    # calculating tendency for theta (first converting back to non-dry theta
    ptr2np(tend_th_ar, size_x, size_z)[1:-1, :] = - (
      ptr2np(th_ar, size_x, size_z)[1:-1, :] -   # old
      arrays["thetad"] # new
#      th_dry2kid(arrays["thetad"], arrays["qv"]) # new
    ) / dt #TODO: check if dt needed


    # calculating tendency for qv
    ptr2np(tend_qv_ar, size_x, size_z)[1:-1, :] = - (
      ptr2np(qv_ar, size_x, size_z)[1:-1, :] - # old                
      arrays["qv"]                             # new 
    ) / dt #TODO: check if dt needed    


    # diagnostics
    if last_diag < it_diag:
      dg.diagnostics(prtcls, arrays, it_diag, size_x, size_z, timestep == 0)
      last_diag = it_diag

    timestep += 1
  except:
    traceback.print_exc()
    return False
  else:
    return True
    
# storing pointers to Python functions
clib.save_ptr(ptrfname, micro_step)

# running Fortran stuff
# note: not using command line arguments, namelist name hardcoded in
#       kid_a_setup/namelists/SC_2D_input.nml 
flib.__main_MOD_main_loop()
