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
Variable_name_l = ["theta", "vapor", "RH", "RH_lib_post_cond", "T_lib_post_cond", "w", "dtheta_mphys", "dqv_mphys", "cloud_number_r20um", "rain_number_r20um", "cloud_mass_r20um", "rain_mass_r20um"]

prsr = ArgumentParser(add_help=True, description='TODO')
prsr.add_argument('--outdir', default="", help='output directory from kid_a_setup/output')
prsr.add_argument('--it_l', nargs='+', type=int, required=None, help='time indexes for plotting, e.g., --it_l 0 2 -1 ')
args = prsr.parse_args()


# reading variables from the netcdf file
def reading_netcdf(netcdf_file, var_l):
    var_d = {}
    for var in var_l + ["z", "x", "time"]:
        var_d[var] = netcdf_file.variables[var][:]
    return var_d


def contour_plot(outdir_path, var_name_l, var_d, it, nr_fig):
    fig = plt.figure(nr_fig, figsize = (8,8))
    fig.suptitle("time = " + str(var_d["time"][it]/60) + " min") 
    x_range = var_d["x"][:-1]
    z_range = var_d["z"][:-1]
    X, Y = np.meshgrid(x_range, z_range)
    nr_pl = 1
    for var in var_name_l:
        print var
        ax = plt.subplot(2,2,nr_pl)
        var_domain = var_d[var][it,:-1,:-1]
        var_min, var_max = var_domain.min(), var_domain.max()
        if var_min == 0.:
            levels_var = np.linspace(var_max * 0.1, var_max, 6)
        else:
            levels_var = np.linspace(var_min, var_max, 6)
        CS = plt.contourf(X, Y, var_domain,  cmap=plt.cm.Blues, alpha=0.7, levels=levels_var)
        plt.xlabel(var + "; min = " +  '%s' % float('%.3g' % var_domain.min()) + 
                   ", max = " + '%s' % float('%.3g' % var_domain.max()), fontsize=10)
        if nr_pl in [1, 3]:
            plt.ylabel(r'height $[m]$', fontsize=10)
        nr_pl += 1

    plt.savefig(os.path.join(outdir_path, "contour_" + str(var_d["time"][it]) + "-".join(var_name_l) + ".pdf"))
    plt.show()



def main(outdir_path, filename="SC_2D_out.nc", variable_name_l=Variable_name_l):
    nf = netcdf.netcdf_file(os.path.join(outdir_path, filename), 'r')
    var_d = reading_netcdf(nf, variable_name_l)
    variable_name_lpl = [variable_name_l[i:i+4] for i in xrange(0, len(variable_name_l), 4)]
    it_l = args.it_l if args.it_l!=None else range(0, var_d["time"].shape[0], 12)
    nr_fig = 1
    for it in it_l:
        for var_name in variable_name_lpl:
            contour_plot(outdir_path, var_name, var_d, it, nr_fig)
            nr_fig += 1

main(os.path.join("../kid_a_setup/output", args.outdir))

