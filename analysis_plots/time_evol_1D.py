from scipy.io import netcdf
import numpy as np
# for skua - TODO
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
from argparse import ArgumentParser
import os
import pdb

# names of variable to plot
Variable_name_l = ["rain_water_path_r20", "rain_water_path_r25", "rain_water_path_r32", \
                   "total_autoconversion_rate_r20", "total_autoconversion_rate_r25","total_autoconversion_rate_r32", \
                   "total_accretion_rate_r20", "total_accretion_rate_r25","total_accretion_rate_r32", \
                   "total_rain_evaporation_rate_r20", "total_rain_evaporation_rate_r25","total_rain_evaporation_rate_r32", \
                   "liquid_water_path", "surface_precipitation_rate"
                  ]
#Variable_name_l = ["vapour", "liquid_drops_number", "liquid_drops_mass", "accumulated_accr20_mass", "accumulated_acnv20_mass", "accumulated_accr32_mass", "accumulated_acnv32_mass", "liquid_water_path", "rain_water_path_r20um"]
#Variable_name_l = ["vapour", "liquid_drops_mass", "liquid_drops_number", "cloud_mass_r20um", "rain_mass_r20um", "accumulated_accr20_mass", "accumulated_acnv20_mass", "accumulated_accr32_mass", "accumulated_acnv32_mass", "liquid_water_path"]
#Variable_name_l = ["vapour", "cloud_mass_r20um", "rain_mass_r20um"]
#Variable_plot_l = Variable_name_l + ["total_water_mass"]#, "cloud+rain_mass", "LWP"]
Variable_plot_l = Variable_name_l# + ["total_water_mass"]#, "cloud+rain_mass", "LWP"]

prsr = ArgumentParser(add_help=True, description='TODO')
prsr.add_argument('--outdir', default="", help='output directory from kid_a_setup/output')
args = prsr.parse_args()

# reading variables from the netcdf file
def reading_netcdf(netcdf_file, var_l):
    var_d = {}
    var_u = {}
    for var in var_l + ["z", "x", "time"]:
        var_d[var] = netcdf_file.variables[var][:]
        var_u[var] = netcdf_file.variables[var].units
    return var_d, var_u


def time_evolution(var_name_l, var_d, var_u, nr_subpl, nr_fig):
    plt.figure(nr_fig, figsize = (20,20))
    nr_pl = 1
    for var in var_name_l:
        print var
        #pdb.set_trace()
        ax = plt.subplot(nr_subpl / 3 + 1,3,nr_pl)
        if(var_d[var].ndim==2):
          var_sum = var_d[var][:,:].sum(axis=1)
        else:
          var_sum = var_d[var][:]
        ax.plot(var_d["time"][1:], var_sum[1:]) #TODO: -999 in the last place
        nr_pl += 1
        plt.ylabel(var + '[' + var_u[var] + ']', fontsize=10)
        plt.xlabel(r'time [s]', fontsize=10)
        
    #plt.savefig("time_evolution_" + "-".join(var_name_l) + ".pdf")
    plt.show()



def main(filename, variable_name_l=Variable_name_l, variable_plot_l=Variable_plot_l, nr_subpl=len(Variable_plot_l)):
    nf = netcdf.netcdf_file(filename, 'r')
    var_d, var_u = reading_netcdf(nf, variable_name_l)
    #var_d["total_water_mass"] = var_d["vapour"] + var_d["cloud_mass_r20um"] + var_d["rain_mass_r20um"]
#    var_d["total_water_mass"] = var_d["vapour"] + var_d["liquid_drops_mass"]
    #var_d["LWP"] = var_d["liquid_drops_mass"] * 25. # dz = 25 m
    #var_d["cloud+rain_mass"] = var_d["cloud_mass_r20um"] + var_d["rain_mass_r20um"]
    variable_plot_lpl = [variable_plot_l[i:i+nr_subpl] for i in xrange(0, len(variable_plot_l), nr_subpl)]
    nr_fig = 1
    for var_name in variable_plot_lpl:
        time_evolution(var_name, var_d, var_u, nr_subpl, nr_fig)
        nr_fig += 1

main(os.path.join(args.outdir, "1D_out.nc"))
