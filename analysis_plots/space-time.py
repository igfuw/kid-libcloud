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
#Variable_name_l = ["theta", "vapour", "RH", "number_of_SDs", "aerosol_number", "w", "dtheta_mphys", "dqv_mphys", "cloud_number_r20um", "rain_number_r20um", "cloud_mass_r20um", "rain_mass_r20um", "liquid_drops_number", "liquid_drops_mean_r", "liquid_drops_mass"]
Variable_name_l = [ \
  "cloud_water_r20", "cloud_water_r25", "cloud_water_r32", \
  "rain_water_r20", "rain_water_r25", "rain_water_r32"  \
  ]

prsr = ArgumentParser(add_help=True, description='TODO')
prsr.add_argument('--outdir', default="", help='output directory from kid_a_setup/output')
prsr.add_argument('--it_l', nargs='+', type=int, required=None, help='time indexes for plotting, e.g., --it_l 0 2 -1 ')
args = prsr.parse_args()


# reading variables from the netcdf file
def reading_netcdf(netcdf_file, var_l):
    var_d = {}
    # comparing with _ante_cond , because post_cond include change in theta due to cond, which is not included in KiD calculations
    # also, for tihs comparison to be correct, save_diag in KiD has to be called before applying mphys and adve dtheta/dqv (i.e. before step_column)
    # but then all the diags in KiD are "pre step", which is not the default behaviour
    var_d["RH_diff"] = -(netcdf_file.variables["RH"][:] - netcdf_file.variables["RH_lib_ante_cond"][:])
    var_d["T_diff"] = -(netcdf_file.variables["temperature"][:] - netcdf_file.variables["T_lib_ante_cond"][:])# / netcdf_file.variables["temperature"][:]
    var_d["p_diff"] = -(netcdf_file.variables["pressure"][:] - netcdf_file.variables["pressure_lib_ante_cond"][:]/1e4)# / netcdf_file.variables["pressure"][:]
    var_d["psat_diff"] = -(netcdf_file.variables["psat"][:] - netcdf_file.variables["psat_lib_formula_ante_cond"][:]/1e2)# / netcdf_file.variables["psat"][:]
    for var in var_l + ["z", "x", "time"]:
        var_d[var] = netcdf_file.variables[var][:]
    return var_d


def contour_plot(outdir_path, var_name_l, var_d, nr_fig):
    fig = plt.figure(nr_fig, figsize = (8,8))
    z_range = var_d["z"][1:]
    x_range = var_d["time"][:]
    print "x_range: ", x_range
    print "z_range: ", z_range
    #pdb.set_trace()
    X, Y = np.meshgrid(x_range, z_range)
    nr_pl = 1
    for var in var_name_l:
        print var
        ax = plt.subplot(2,2,nr_pl)
        var_domain = var_d[var][:,1:]
        var_min, var_max = var_domain.min(), var_domain.max()
        if var_min < var_max:
            if var_min == 0.: 
                levels_var = np.linspace(var_max * 0.1, var_max, 20)
            else:
                levels_var = np.linspace(var_min, var_max, 20)
            print levels_var
            CS = plt.contourf(X, Y, var_domain.transpose(),  cmap=plt.cm.Blues, alpha=0.7, levels=levels_var)
            plt.xlabel(var + "; min = " +  '%s' % float('%.3g' % var_domain.min()) +
                       ", max = " + '%s' % float('%.3g' % var_domain.max()) + "\n" +
                       "min_level = " +  '%s' % float('%.3g' % levels_var[0]),
                       fontsize=10)
            if nr_pl in [1, 3]:
                plt.ylabel(r'height $[m]$', fontsize=10)
        nr_pl += 1

    plt.savefig(os.path.join(outdir_path, "contour_".join(var_name_l) + ".pdf"))
    plt.show()



def main(outdir_path, filename="1D_out.nc", variable_name_l=Variable_name_l):
    nf = netcdf.netcdf_file(os.path.join(outdir_path, filename), 'r')
    var_d = reading_netcdf(nf, variable_name_l)
    variable_name_l  += ["RH_diff", "T_diff", "p_diff", "psat_diff"]
    variable_name_lpl = [variable_name_l[i:i+4] for i in xrange(0, len(variable_name_l), 4)]
    nr_fig = 1
    for var_name in variable_name_lpl:
        contour_plot(outdir_path, var_name, var_d, nr_fig)
        nr_fig += 1

main(args.outdir)
