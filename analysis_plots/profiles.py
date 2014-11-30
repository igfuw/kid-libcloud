from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker

# reading variables from the netcdf file
def reading_netcdf(file_r, var):
    nf = netcdf.netcdf_file(file_r, 'r')
    var_nctab = nf.variables[var][:]
    return var_nctab
