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
Variable_name_l = ["theta", "vapor", "RH", "RH_lib_post_cond", "aerosol_number", "w", "dtheta_mphys", "dqv_mphys",  "cloud_number_r20um", "rain_number_r20um", "cloud_mass_r20um", "rain_mass_r20um"]

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


def plotting_profiles(outdir_path, var_name_l, var_d, it_plot_l):
    plt.figure(1, figsize = (8,8))
    nr_pl = 1
    for var in var_name_l:
        print var
        ax = plt.subplot(2,2,nr_pl)
        legend = []
        for it in it_plot_l:
            legend.append("time = " + str(var_d["time"][it]/60) + " min")
            var_hor_av = np.mean(var_d[var][it,:,:], axis=1)
            print "min, max", var_hor_av.min(), var_hor_av.max()
            ax.plot(var_hor_av[:-1], var_d["z"][:-1]) #TODO: -999 in the last place, using missing arrays
        nr_pl += 1
        #ax.set_xlim(var_d[var][1:,:-1,:].min(), var_d[var][1:,:-1:,:].max())
        print len(ax.xaxis.get_major_ticks())
        for i, tick in enumerate(ax.xaxis.get_major_ticks()):
            if i % 2 != 0:
                tick.label1On = False
        for item in plt.xticks()[1] + plt.yticks()[1]:
            item.set_fontsize(7)
        plt.xlabel(var, fontsize=10)
        plt.ylabel(r'height $[m]$', fontsize=10)
        if nr_pl == 3:
            plt.legend(legend, prop = FontProperties(size=6))

    plt.savefig(os.path.join(outdir_path, "profiles_" + "-".join(var_name_l) + ".pdf"))
    plt.show()



def main(outdir_path, filename="SC_2D_out.nc", variable_name_l=Variable_name_l):
    nf = netcdf.netcdf_file(os.path.join(outdir_path, filename), 'r')
    var_d = reading_netcdf(nf, variable_name_l)
    variable_name_lpl = [variable_name_l[i:i+4] for i in xrange(0, len(variable_name_l), 4)]
    it_l = args.it_l if args.it_l!=None else range(0, var_d["time"].shape[0], 4)
    for var_name in variable_name_lpl:
        plotting_profiles(outdir_path, var_name, var_d, it_plot_l=it_l)


main(os.path.join("../kid_a_setup/output", args.outdir))
