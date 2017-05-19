__author__ = 'david'

import numpy as np
import netCDF4
import pylab
import os
import shutil
import copy
import glob
import subprocess
from scipy.io import netcdf as nc

#
# This module contains some useful utilities:
#
# 1) log_write: write to the log file
#
# 2) read_nml: read in the namelist for the optimization
#
# 3) run_forward: integrates the forward model for each optimization timestep
#

###################################
#                                 #
# Write a message to the log file #
#                                 #
###################################

def log_write(path_output, string, create = 'false'):
	
	if create == 'true':
		log_full = open(path_output + 'LOGFILE.txt','w')
	else:
		log_full = open(path_output + 'LOGFILE.txt','a')

	log_full.write(string)
	log_full.write('\n')
	log_full.close()




def test(a):
	print a


#################
#               #
# Read namelist #
#               #
#################

def read_nml(path_scripts):

	g = open(path_scripts + 'namelist.CG.txt','r')
	lines = g.readlines()
	for line in lines:
		if (len(line) > 0):
			a = line.split()[0]
			b = line.split()[2]
			if (a == 'INPUT_DIR'):
				INPUT_DIR = b
			elif (a == 'FORWARD_RUN_DIR'):
				FORWARD_RUN_DIR = b
			elif (a == 'ADJ_RUN_DIR'):
				ADJ_RUN_DIR = b
			elif (a == 'OUTPUT_DIR'):
				OUTPUT_DIR = b
			elif (a == 'Nt'):
				Nt = int(b)
			elif (a == 'dx'):
				dx = float(b)
			elif (a == 'dy'):
				dy = float(b)
			elif (a == 'P_scale'):
				P_scale = float(b)
			elif (a == 'U_scale'):
				U_scale = float(b)
			elif (a == 'V_scale'):
				V_scale = float(b)
			elif (a == 'W_scale'):
				W_scale = float(b)
			elif (a == 'MU_scale'):
				MU_scale = float(b)
			elif (a == 'PH_scale'):
				PH_scale = float(b)
			elif (a == 'QVAPOR_scale'):
				QVAPOR_scale = float(b)
			elif (a == 'max_iter'):
				max_iter = int(b)
			elif (a == 'tol'):
				tol = float(b)
			elif (a == 'storm_name'):
				storm_name = b
			elif (a == 'opt_interval'):
				opt_interval = int(b)
			elif (a == 'start_year'):
				start_year = int(b)
			elif (a == 'start_day'):
				start_day = int(b)
			elif (a == 'start_month'):
				start_month = int(b)
			elif (a == 'start_hour'):
				start_hour = int(b)
			elif (a == 'phys_test'):
				phys_test = str(b)
			elif (a == 'MET_DIR'):
				MET_DIR = str(b)
			elif (a == 'run_storm'):
				run_storm = str(b)
			elif (a == 'alpha_min'):
				alpha_min = float(b)
			elif (a == 'R'):
				R = float(b)
			elif (a == 'end'):
				break

	scales = [P_scale, U_scale, V_scale, W_scale, PH_scale, MU_scale, QVAPOR_scale]
	dirs = [path_scripts, INPUT_DIR, FORWARD_RUN_DIR, ADJ_RUN_DIR, OUTPUT_DIR, MET_DIR]
	params = [Nt, dx, dy, alpha_min, R]
	starts = [start_year, start_month, start_day, start_hour]
	g.close()
	
	return scales, dirs, params, starts, max_iter, tol, storm_name, opt_interval, phys_test, run_storm



def run_forward(Nt, rsts, dirs, count):
	
	path_input = dirs[1]
	path_forward = dirs[2]
	path_adj = dirs[3]
	path_output = dirs[4]

        log_write(path_output, 'Calling forward model')

	for i in range(Nt):

		if i == 0:

			os.chdir(path_forward)
			shutil.copy2(path_input + 'namelist.input.forward_0', path_forward + 'namelist.input')
			os.system('mpirun -np 28 ./wrf.exe')
			shutil.move(path_forward + 'rsl.out.0000', path_output + 'rsl.' + str(count) + '.' + str(i))
			for filename in glob.glob(os.path.join(path_forward, 'wrfout*')):
				shutil.move(filename, path_output + 'wrfout_' + str(i) + '.nc')

		else:

			os.chdir(path_forward)
			shutil.copy2(path_input + rsts[i-1] +'.nc', path_forward + rsts[i-1])
			shutil.copy2(path_input + 'namelist.input.forward_' + str(i), path_forward + 'namelist.input')
			os.system('mpirun -np 28 ./wrf.exe')
			os.remove(rsts[i-1])
			shutil.move(path_forward + 'rsl.out.0000', path_output + 'rsl.' + str(count) + '.' + str(i))
			for filename in glob.glob(os.path.join(path_forward, 'wrfout*')):
				shutil.move(filename, path_output + 'wrfout_' + str(i) + '.nc')
	


