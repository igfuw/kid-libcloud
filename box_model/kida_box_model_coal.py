# test coalescence algorithm for geometric kernel with Hall efficiencies for drops with r>30um and Davis&Rogers efficiencies for smaller ones
# by comparing mass density function with results of EFM modeling

import sys
import time

sys.path.insert(0, "/home/piotr/usr/local/lib/python3/dist-packages/")

from libcloudphxx import lgrngn
from math import exp, sqrt
import numpy as np

from scipy.stats import gamma
import matplotlib.pyplot as plt

import netCDF4

#fig, axarr = plt.subplots(1, 2)

#axarr[0].set(xscale='log', yscale='log', xlabel='diameter[um]', ylabel='number density function over ln(d) [1/cc]')
#axarr[1].set(xscale='log', yscale='log', xlabel='diameter[um]', ylabel='mass density function over ln(d) [g/kg]')

# initial values for the search of shape parameter
init_theta = {2.5: {50e6: 3.8e-6, 150e6 : 2.6e-6, 300e6: 2.1e-6},
              0:   {50e6: 9.3e-6, 150e6 : 6.462821053522438e-06, 300e6: 5.129933456758218e-06}}


#total time of simulation
simulation_time = 3600
output_interval = 2 # dsd stored each 2s

# total number of SD in the whole ensemble
tot_sd_no = 1e5

#root mean square deviation
def RMSD(a1, a2):
  nonempty = 0 
  tot = 0 
  for i in range(a1.size):
    if(a1[i] > 0 or a2[i] > 0): 
      tot+=pow(a1[i] - a2[i], 2)
      nonempty+=1
  return np.sqrt(tot/nonempty)

# numerical integral of the distirbution to get total mass
def int_gamma(alfa, theta):
  r_edges = np.linspace(gamma(alfa, scale = theta).ppf(0.0001),
                  gamma(alfa, scale = theta).ppf(0.9999), 10001) # bin edges
                  
#  print(r_edges)

  dr = np.zeros(len(r_edges)-1)
  r_centers = np.zeros(len(r_edges)-1)
  mass = np.zeros(len(r_edges)-1)
  conc = np.zeros(len(r_edges)-1)
  
  dr[:] = r_edges[1:] - r_edges[:-1]
  r_centers[:] = r_edges[:-1] + 0.5 * dr[:]
  water_dens = 1e3 # kg/m^3
  conc = n_zero * dr * gamma(alfa, scale = theta).pdf(r_centers)
  mass = conc * 4./3. * np.pi * r_centers**3 * water_dens
  return np.sum(conc), np.sum(mass)

opts_init = lgrngn.opts_init_t()
opts_init.dt = output_interval


opts_init.dx = 1
opts_init.dz = 1
opts_init.nz = 1
opts_init.z1 = opts_init.dz * opts_init.nz

kappa = 1e-10


opts_init.kernel = lgrngn.kernel_t.hall
#opts_init.kernel = lgrngn.kernel_t.hall_davis_no_waals
opts_init.terminal_velocity = lgrngn.vt_t.beard76

opts_init.sd_conc_large_tail = 1
opts_init.aerosol_independent_of_rhod = 1

opts = lgrngn.opts_t()
opts.adve = False
opts.sedi = False
opts.cond = False
opts.coal = True
opts.chem = False

# calc output bin edges, in diameter, mass doubling
d0 = 3.125e-6 # first bin left edge
n_bins = 50
bins_m = np.zeros(n_bins+1)
bins_m[0] = d0**3 # mass of droplets with d0, without 1/6 pi liq_dens
for i in np.arange(1,n_bins+1):
  bins_m[i] = 2*bins_m[i-1]
bins_d = pow(bins_m, 1./3.)
bins_d_ctr = np.zeros(n_bins)
bins_d_wdt = np.zeros(n_bins)
for i in range(n_bins):
  bins_d_ctr[i] = (bins_d[i] + bins_d[i+1]) / 2.
  bins_d_wdt[i] = (bins_d[i+1] - bins_d[i])

res_conc = np.zeros(n_bins)
res_mass = np.zeros(n_bins)

def diag_dsd(conc, mass):
  for i in range(conc.size) :
    prtcls.diag_wet_rng(bins_d[i] / 2., bins_d[i+1] / 2.)
    prtcls.diag_wet_mom(0)
    conc[i]= np.frombuffer(prtcls.outbuf()).mean()

    prtcls.diag_wet_rng(bins_d[i] / 2., bins_d[i+1] / 2.)
    prtcls.diag_wet_mom(3)
    mass[i]= np.frombuffer(prtcls.outbuf()).mean() * 4. / 3. * np.pi * 1e3

def diag_tot_conc_mass():
  prtcls.diag_all()
  prtcls.diag_wet_mom(0)
  conc = np.frombuffer(prtcls.outbuf()).mean()
  prtcls.diag_all()
  prtcls.diag_wet_mom(3)
  mass = np.frombuffer(prtcls.outbuf()).mean() * 4./3. * np.pi * 1e3
  return conc, mass


liq_mass = 0.001 # kg/kg

for n_zero in [50e6, 150e6, 300e6]:
  for shape_param in [0]:
  #for shape_param in [0, 2.5]:
    # initial gamma distribution of droplet radii, defined as in Pruppacher & Klett Eq. (11-108)
    alfa = shape_param + 1 # shape parameter in the shape-rate representation (see https://en.wikipedia.org/wiki/Gamma_distribution)
    #beta = shape_param + 1 # rate parameter in the shape-rate representation (see https://en.wikipedia.org/wiki/Gamma_distribution)
    #theta = 1 / beta       # scale parameter in the shape-scale representation (see https://en.wikipedia.org/wiki/Gamma_distribution)
    theta = init_theta[shape_param][n_zero]
    
    # look for a scale parameter that would give desrired liquid mass
    eps = 1e-3 # relative tolerance
    while True:
      conc, mass = int_gamma(alfa, theta)
      if abs(mass - liq_mass) / liq_mass < eps:
        break
      if mass > liq_mass:
        theta /= 1.0001
      if mass < liq_mass:
        theta *= 1.0001
    
    print(n_zero, shape_param, ': ', theta, conc, mass)
    
    # plotting the gamma distribution, as a function of diameter
    #x = np.linspace(gamma(alfa, scale = theta).ppf(0.01),
    #                gamma(alfa, scale = theta).ppf(0.99), 100)
    #axarr[0].plot(2*x, n_zero * gamma(alfa, scale = theta).pdf(x) / 2., # *2 and /2 because its in radius
    #       'r-', lw=5, alpha=0.6, label='gamma pdf')
    
    # gamma distribution used to init dsd in libcloud. function of ln(r)
    def gammalnr(lnr):
      r=np.exp(lnr)
      if(shape_param==0 and r<1e-10): # shape_param==0 gives non-zero n for r=0, we cut it off
        return 0
      return n_zero * gamma(alfa, scale = theta).pdf(r) * r
    
    opts_init.dry_distros = {kappa:gammalnr}
#    if(shape_param == 0):
#      opts_init.rd_min = 1e-9 # don't init smaller than 0.1nm. For shape_param=0, init distro is nonzero for r=0
    
    
    for sd_conc in [100]:
      for dt_coal in [1]:
        # netcdf output stuff
        ncfile = netCDF4.Dataset('UWLCM_box_N'+str(n_zero)+'_shape'+str(shape_param)+'_SD'+str(sd_conc)+'_dtCoal'+str(dt_coal)+'_totSD'+str(tot_sd_no)+'_dsd.nc',mode='w',format='NETCDF4_CLASSIC') 
        
        #dimensions
        bin_edges_dim = ncfile.createDimension('bin_edges', n_bins+1)
        bin_data_dim = ncfile.createDimension('bin_data', n_bins)
        time_dim = ncfile.createDimension('time', (simulation_time / 2) + 1) # /2 because we output each 2s
        
        #variables
        times = ncfile.createVariable('output_time', np.float32, ('time',))
        times.units = "s"
        times.long_name = "time at which dsd is diagnosed"
        times[:] = np.arange(0,simulation_time+1,2)
        
        bin_edges = ncfile.createVariable('bin_edges', np.float32, ('bin_edges',))
        bin_edges.units = "m"
        bin_edges.long_name = "droplet diameter at the edges of output bins"
        bin_edges[:] = bins_d[:]
        
        bin_conc = ncfile.createVariable('bin_conc', np.float32, ('time','bin_data'))
        bin_conc.units = "m^{-3}"
        bin_conc.long_name = "concentration of droplets with sizes within a bin"
        
        bin_mass = ncfile.createVariable('bin_mass', np.float32, ('time','bin_data'))
        bin_mass.units = "kg m^{-3}"
        bin_mass.long_name = "mass of droplets with sizes within a bin"
    
        opts_init.nx = int(tot_sd_no / sd_conc)
        opts_init.x1 = opts_init.dx * opts_init.nx
    
        rhod =   1. * np.ones((opts_init.nx, opts_init.nz))
        th   = 300. * np.ones((opts_init.nx, opts_init.nz))
        rv   = 0.01 * np.ones((opts_init.nx, opts_init.nz))
    
        opts_init.rng_seed = int(time.time())
        opts_init.sd_conc = sd_conc
        opts_init.n_sd_max = int(opts_init.nx * opts_init.sd_conc * 1.1) # make room for large tail
        opts_init.sstp_coal = int(output_interval / dt_coal)
    
        print('sd_conc: ', sd_conc, 'dt_coal: ', dt_coal)
      
        try:
          prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
        except:
          prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)
        
        prtcls.init(th, rv, rhod)
      
        init_conc, init_mass = diag_tot_conc_mass()
        print("initial conc of droplets: ", init_conc, '[1/kg]')
        print("initial mass of droplets: ", init_mass, '[kg/kg]')
      
        diag_dsd(res_conc, res_mass)
        # res_conc is N, the number of droplets in the disred range of radii. we want to plot number density function, n, over ln(D). by definition N = n(ln(D)) d ln(D), 
        # so n(ln(D)) = N / d ln(D) = N * dD / dD / d ln(D) = N / dD * D
#        axarr[0].plot(bins_d_ctr*1e6, res_conc / bins_d_wdt * bins_d_ctr / 1e6, label=None) # / 1e6 to get 1/cc 
#        axarr[1].plot(bins_d_ctr*1e6, res_mass / bins_d_wdt * bins_d_ctr * 1e3, label=None) # * 1e3 to get grams
    
        # output initial dsd to netcdf
        bin_conc[0,:] = res_conc[:]
        bin_mass[0,:] = res_mass[:]
        
        #simulation loop
        for t in range(int(simulation_time / output_interval)+1):
          print(t*output_interval / simulation_time * 100, ' %', end='\r')
          prtcls.step_sync(opts, th, rv, rhod)
          prtcls.step_async(opts)
          # output dsd to netcdf
          diag_dsd(res_conc, res_mass)
          bin_conc[t,:] = res_conc[:]
          bin_mass[t,:] = res_mass[:]
            
      #  axarr[0].plot(bins_d_ctr*1e6, res_conc / bins_d_wdt * bins_d_ctr / 1e6, label='sd_conc: '+str(sd_conc)+' dt_coal: '+str(dt_coal)) 
      #  axarr[1].plot(bins_d_ctr*1e6, res_mass / bins_d_wdt * bins_d_ctr * 1e3, label='sd_conc: '+str(sd_conc)+' dt_coal: '+str(dt_coal)) 
      
        final_conc, final_mass = diag_tot_conc_mass()
        print("final conc of droplets: ", final_conc, '[1/kg]')
        print("final mass of droplets: ", final_mass, '[kg/kg]')
        print("final mass of droplets in outpput bins: ", np.sum(res_mass), '[kg/kg]')
      
 #   axarr[0].legend()
 #   axarr[1].legend()
 #   plt.show()

    ncfile.close()
