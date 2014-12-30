from scipy.io import netcdf
import numpy as np
# for skua - TODO
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import pdb

# names of variable to plot
Variable_name_l = ["theta", "vapor", "RH", "RH_lib_post_cond", "T_lib_post_cond", "w", "dtheta_mphys", "dqv_mphys",  "cloud_number_r20um", "rain_number_r20um", "cloud_mass_r20um", "rain_mass_r20um"]

# reading variables from the netcdf file
def reading_netcdf(netcdf_file, var_l):
    var_d = {}
    for var in var_l + ["z", "x", "time"]:
        var_d[var] = netcdf_file.variables[var][:]
    return var_d


def plotting_profiles(var_name_l, var_d):
    plt.figure(1, figsize = (8,8))
    nr_pl = 1
    for var in var_name_l:
        print var
        ax = plt.subplot(2,2,nr_pl)
        legend = []
        for it in range(var_d[var].shape[2]):
            legend.append("time = " + str(var_d["time"][it]))
            var_hor_av = np.mean(var_d[var][:,:,it], axis=1)
            ax.plot(var_hor_av[:-1], var_d["z"][:-1]) #TODO: -999 in the last place
        nr_pl += 1
        plt.xlabel(var, fontsize=10)
        plt.ylabel(r'height $[m]$', fontsize=10)
        plt.legend(legend, prop = FontProperties(size=8))

    plt.savefig("profiles_" + "-".join(var_name_l) + ".pdf")
    plt.show()



def main(filename, variable_name_l=Variable_name_l):
    nf = netcdf.netcdf_file(filename, 'r')
    var_d = reading_netcdf(nf, variable_name_l)
    variable_name_lpl = [variable_name_l[i:i+4] for i in xrange(0, len(variable_name_l), 4)]
    for var_name in variable_name_lpl:
        plotting_profiles(var_name, var_d)


main("../kid_a_setup/output/SC_2D_out.nc")
