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
Variable_name_l = ["theta", "vapor", "RH", "RH_lib_post_cond", "T_lib_post_cond", "w", "dtheta_mphys", "dqv_mphys", "cloud_number_r20um", "rain_number_r20um", "cloud_mass_r20um", "rain_mass_r20um"]

# reading variables from the netcdf file
def reading_netcdf(netcdf_file, var_l):
    var_d = {}
    for var in var_l + ["z", "x", "time"]:
        var_d[var] = netcdf_file.variables[var][:]
    return var_d


def contour_plot(var_name_l, var_d, it):
    plt.figure(1, figsize = (8,8))
    x_range = var_d["x"][:-1]
    z_range = var_d["z"][:-1]
    X, Y = np.meshgrid(x_range, z_range)
    nr_pl = 1
    for var in var_name_l:
        print var
        ax = plt.subplot(2,2,nr_pl)
        #legend = []
        #legend.append("time = " + str(var_d["time"][it]))
        var_domain = var_d[var][:-1,:-1,it]
        CS = plt.contourf(X, Y, var_domain,  cmap=plt.cm.Blues, alpha=0.7)
        nr_pl += 1
        plt.xlabel(var + "; min = " +  '%s' % float('%.3g' % var_domain.min()) + 
                   ", max = " + '%s' % float('%.3g' % var_domain.max()), fontsize=10)
        plt.ylabel(r'height $[m]$', fontsize=10)
#plt.legend(legend, prop = FontProperties(size=8))

    plt.savefig("contour_" + str(var_d["time"][it]) + "-".join(var_name_l) + ".pdf")
    plt.show()



def main(filename, variable_name_l=Variable_name_l, it_l=[0,-1]):
    nf = netcdf.netcdf_file(filename, 'r')
    var_d = reading_netcdf(nf, variable_name_l)
    variable_name_lpl = [variable_name_l[i:i+4] for i in xrange(0, len(variable_name_l), 4)]
    for it in it_l:
        for var_name in variable_name_lpl:
            contour_plot(var_name, var_d, it)


main("../kid_a_setup/output/SC_2D_out.nc")
