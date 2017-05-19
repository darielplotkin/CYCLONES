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

#
# This module contains the routine the iterates the optimization:
#
# 1) iterate: iterates the minimization given a tolerance for
#             the difference in cost function between steps
#

###########################################
#                                         #
# Iterate the optimization to convergence #
#                                         #
###########################################

def iterate(restarts, dirs, params, scales, tol, max_iter, storm_name):

	path_input = dirs[1]
	path_forward = dirs[2]
	path_adj = dirs[3]
	path_output = dirs[4]
	path_met = dirs[5]
	Nt = params[0]
	alpha_min = params[3]
	R = params[4]

	count = 0
	dif = 2*tol

	while (dif > tol and count < max_iter):
	
		alpha = 1.
		
		if count == 0:

			utils.run_forward(Nt, restarts, dirs, count)

			U_list = variables.x_to_U(restarts, Nt, scales, dirs)
			g_old = cost.get_g(Nt, R, scales, U_list, dirs, restarts, storm_name)
			g_new = 2*g_old

		if count > 0:

			utils.log_write(path_output, '#########################################')
			utils.log_write(path_output, 'Beginning optimization iteration number '+str(count))
			utils.log_write(path_output, '#########################################')	

			if (count == 1):
				[p_old_list, grad_old_list, U_list, g_new] = update.update_U(params, scales, restarts, storm_name, U_list, dirs, count)

			else:
				[p_old_list, grad_old_list, U_list, g_new] = update.update_U(params, scales, restarts, storm_name, U_list, dirs, count, p_old_list, grad_old_list)

				utils.log_write(path_output, 'counter is: ' +str(count) + '     value function is: ' + str(g_old))

		dif = np.abs(g_old - g_new)
		g_old = g_new

		count += 1
