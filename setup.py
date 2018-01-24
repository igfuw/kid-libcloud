import numpy as np
import libcloudphxx as libcl


# dictionary of simulation parameters
params = {
  "real_t" : np.float64,
  "sd_conc" : 16,
  "spinup_rain" : 60*30, # number of timesteps during which sedimentation and coalescence are off
  "spinup_smax" : 60*10,
  "kappa" : .61,
  "meanr" : .04e-6,
  "gstdv" : 1.4,
  "n_tot" : 100e6,
  "n_bins": 34,              # \__ from the TAU example file @ KiD-A website                   
  "bin0_D_upper" : 3.125e-6, # /                                                               
  # bins considered as cloud water (all below is aerosol, all above is rain) 
  "bins_qc_r20um" : np.arange(1, 12), # r=20 um threshold setting (actually D=39.7 um)
  "bins_qc_r32um" : np.arange(1, 14), # r=32 um threshold setting (actually D=63.0 um)
  "sstp_cond" : 10,
  "sstp_coal" : 10
}

opts = libcl.lgrngn.opts_t()
opts.cond = True
opts.coal = True
opts.adve = True
opts.sedi = True
opts.chem = False
opts.RH_max = 1.005
opts.kernel = libcl.lgrngn.kernel_t.geometric
