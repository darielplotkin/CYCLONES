__author__ = 'david'

import numpy as np
import netCDF4
import pylab
import os
import shutil
import copy
import glob
import subprocess
import nc_copy
from scipy.io import netcdf as nc
import init
import variables
import update
import cost
import utils
import finalize
import iterate




###############
#             #
# Main Script #
#             #
###############

##### get current (scripts) directory ########
path_scripts = os.getcwd()+'/'

###### read the namelist to get parameters ###########
scales, dirs, params, starts, max_iter, tol, storm_name, opt_interval, phys_test, run_storm = utils.read_nml(path_scripts)

######### Initialize domains and al necessary files for optimization #############
restarts = init.initialize(dirs, params, storm_name, opt_interval, starts, phys_test, run_storm)

######## Iterate the optimization ########
iterate.iterate(restarts, dirs, params, scales, tol, max_iter, storm_name)

############ create a netCDF file with the final trajectory ##############
finalize.post_process(params, dirs, restarts, opt_interval)

