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
import utils

#
# This module contains routines for creating 
# all of the files used in action minimization.
# Routines are:
#
# 1) initial_profile: creates an interpolated 
#    initial trajectory given start and end states
#
# 2) physical_profile: creates a physical initial
#    trajectory given start state
#
# 3) physical_profile: creates a physical initial
#    trajectory with a small perturbation given start state
#
# 4) get_namelists: creates all necessary wrf namelists
#    for forward and adjoint integration
#
# 5) final_sens_init: creates final sensitivity files
#    for adjoint integration of WRF model
#
# 6) initialize: runs all of the above
#

################################
#                              #
# Create an initial trajectory #
# interpolated between initial #
# and final states             #
#                              #
################################


def initial_profile(Nt, path_forward, path_input, path_output, storm_name, opt_interval, starts, run_storm):

	if run_storm == 'true':	
		shutil.copy2(path_input + 'namelist.input.forward_long', path_forward + 'namelist.input')
		os.chdir(path_forward)
		os.system('mpirun -np 28 ./wrf.exe')

		for filename in glob.glob(os.path.join(path_input, 'wrfout_d01*')):
			shutil.move(filename,  path_input + 'wrfout_' + storm_name + '.nc')
	
	else:
		shutil.copy2(path_input + 'namelist.input.forward_long_abbr', path_forward + 'namelist.input')
		os.chdir(path_forward)
		os.system('mpirun -np 28 ./wrf.exe')

	restarts = []
	
	for i in range(Nt+1):

		if i == 0:
			d1 = starts[2]/10
			d2 = starts[2]%10

			h1 = starts[3]/10
			h2 = starts[3]%10			
		else:
			h2 += opt_interval
			if h2 >= 10:
				h1 += h2/10
				h2 = h2%10
			if (h1 == 2 and h2 == 4):
				h1 = 0
				h2 = 0
				d2 += 1
				if d2%10 == 0:
					d1 += 1 
					d2 = 0
			
			cur_day = str(d1) + str(d2)
			cur_hour = str(h1) + str(h2)

			if starts[1] < 10:
				start_month = '0' + str(starts[1])
			else:
				start_month = str(starts[1])
				
			string = 'wrfrst_d01_' + str(starts[0]) + '-' + start_month + '-' + cur_day + '_' + cur_hour + ':00:00'
			restarts.append(string)
			shutil.move(path_input + string, path_input + string + '.nc')
			shutil.copy2(path_input + string + '.nc', path_input + 'ORIG_' + string + '.nc')

	files = []
	file_end = netCDF4.Dataset(path_input + 'wrfout_' + storm_name + '.nc','r+')
	file_init = netCDF4.Dataset(path_input + 'wrfinput_d01_IC.nc')

	SST_final = file_end['SST'][70]
	CANWAT_final = file_end['CANWAT'][70]
	MUB_final = file_end['MUB'][70]
	P_HYD_final = file_end['P_HYD'][70]
	PB_final = file_end['PB'][70]
	PHB_final = file_end['PHB'][70]
	PSFC_final = file_end['PSFC'][70]
	Q2_final = file_end['Q2'][70]
	QCLOUD_final = file_end['QCLOUD'][70]
	QRAIN_final = file_end['QRAIN'][70]
	T_final = file_end['T'][70]
	T2_final = file_end['T2'][70]
	TH2_final = file_end['TH2'][70]
	TSK_final = file_end['TSK'][70]
	U10_final = file_end['U10'][70]
	V10_final = file_end['V10'][70]
	P_final = file_end['P'][70]
	U_final = file_end['U'][70]
	V_final = file_end['V'][70]
	W_final = file_end['W'][70]
	MU_final = file_end['MU'][70]
	PH_final = file_end['PH'][70]
	QVAPOR_final = file_end['QVAPOR'][70]

	SST_init = file_init['SST'][0]
	CANWAT_init = file_init['CANWAT'][0]
	MUB_init = file_init['MUB'][0]
	P_HYD_init = file_init['P_HYD'][0]
	PB_init = file_init['PB'][0]
	PHB_init = file_init['PHB'][0]
	PSFC_init = file_init['PSFC'][0]
	Q2_init = file_init['Q2'][0]
	QCLOUD_init = file_init['QCLOUD'][0]
	QRAIN_init = file_init['QRAIN'][0]
	T_init = file_init['T'][0]
	T2_init = file_init['T2'][0]
	TH2_init = file_init['TH2'][0]
	TSK_init = file_init['TSK'][0]
	U10_init = file_init['U10'][0]
	V10_init = file_init['V10'][0]
	P_init = file_init['P'][0]
	U_init = file_init['U'][0]
	V_init = file_init['V'][0]
	W_init = file_init['W'][0]
	MU_init = file_init['MU'][0]
	PH_init = file_init['PH'][0]
	QVAPOR_init = file_init['QVAPOR'][0]

	SST_interval = (SST_final - SST_init)/Nt
	CANWAT_interval = (CANWAT_final - CANWAT_init)/Nt
	MUB_interval = (MUB_final - MUB_init)/Nt
	P_HYD_interval = (P_HYD_final - P_HYD_init)/Nt
	PB_interval = (PB_final - PB_init)/Nt
	PHB_interval = (PHB_final - PHB_init)/Nt
	PSFC_interval = (PSFC_final - PSFC_init)/Nt
	Q2_interval = (Q2_final - Q2_init)/Nt
	QCLOUD_interval = (QCLOUD_final - QCLOUD_init)/Nt
	QRAIN_interval = (QRAIN_final - QRAIN_init)/Nt
	T_interval = (T_final - T_init)/Nt
	T2_interval = (T2_final - T2_init)/Nt
	TH2_interval = (TH2_final - TH2_init)/Nt
	TSK_interval = (TSK_final - TSK_init)/Nt
	U10_interval = (U10_final - U10_init)/Nt
	V10_interval = (V10_final - V10_init)/Nt
	P_interval = (P_final - P_init)/Nt
	U_interval = (U_final - U_init)/Nt
	V_interval = (V_final - V_init)/Nt
	W_interval = (W_final - W_init)/Nt
	MU_interval = (MU_final - MU_init)/Nt
	PH_interval = (PH_final - PH_init)/Nt
	QVAPOR_interval = (QVAPOR_final - QVAPOR_init)/Nt	
	
	file_end.close()

	for i in range(Nt):

		files.append(netCDF4.Dataset(path_input + restarts[i] + '.nc','r+'))
		t = files[i]['Times'][:]

		if i == 0:
			d1 = starts[2]/10
			d2 = starts[2]%10

			h1 = starts[3]/10
			h2 = starts[3]%10			
		else:
			h2 += opt_interval
			if h2 >= 10:
				h1 += h2/10
				h2 = h2%10
			if (h1 == 2 and h2 == 4):
				h1 = 0
				h2 = 0 
				d2 += 1
				if d2%10 == 0:
					d1 += 1 
					d2 = 0

		t[0][8] = str(d1)
		t[0][9] = str(d2)
		t[0][11] = str(h1)
		t[0][12] = str(h2)

		files[i]['T_INIT'][0] = file_init['T_INIT'][0]
		files[i]['VOCE'][0] = file_init['VOCE'][0]
		files[i]['UOCE'][0] = file_init['UOCE'][0]
		files[i]['SST'][0] = file_init['SST'][0] + SST_interval*(i+1)
		files[i]['CANWAT'][0] = file_init['CANWAT'][0] + CANWAT_interval*(i+1)
		files[i]['MUB'][0] = file_init['MUB'][0] + MUB_interval*(i+1)
		files[i]['P_HYD'][0] = file_init['P_HYD'][0] + P_HYD_interval*(i+1)
		files[i]['PB'][0] = file_init['PB'][0] + PB_interval*(i+1)
		files[i]['PHB'][0] = file_init['PHB'][0] + PHB_interval*(i+1)
		files[i]['PSFC'][0] = file_init['PSFC'][0] + PSFC_interval*(i+1)
		files[i]['Q2'][0] = file_init['Q2'][0] + Q2_interval*(i+1)
		files[i]['QCLOUD'][0] = file_init['QCLOUD'][0] + QCLOUD_interval*(i+1)
		files[i]['QRAIN'][0] = file_init['QRAIN'][0] + QRAIN_interval*(i+1)
		files[i]['T_1'][0] = file_init['T'][0] + T_interval*(i+1)
		files[i]['T_2'][0] = file_init['T'][0] + T_interval*(i+1)
		files[i]['T2'][0] = file_init['T2'][0] + T2_interval*(i+1)
		files[i]['TH2'][0] = file_init['TH2'][0] + TH2_interval*(i+1)
		files[i]['TSK'][0] = file_init['TSK'][0] + TSK_interval*(i+1)
		files[i]['U10'][0] = file_init['U10'][0] + U10_interval*(i+1)
		files[i]['V10'][0] = file_init['V10'][0] + V10_interval*(i+1)
		files[i]['P'][0] = file_init['P'][0] + P_interval*(i+1)
		files[i]['U_1'][0] = file_init['U'][0] + U_interval*(i+1)
		files[i]['U_2'][0] = file_init['U'][0] + U_interval*(i+1)
		files[i]['V_1'][0] = file_init['V'][0] + V_interval*(i+1)
		files[i]['V_2'][0] = file_init['V'][0] + V_interval*(i+1)
		files[i]['W_1'][0] = file_init['W'][0] + W_interval*(i+1)
		files[i]['W_2'][0] = file_init['W'][0] + W_interval*(i+1)
		files[i]['MU_1'][0] = file_init['MU'][0] + MU_interval*(i+1)
		files[i]['MU_2'][0] = file_init['MU'][0] + MU_interval*(i+1)
		files[i]['PH_1'][0] = file_init['PH'][0] + PH_interval*(i+1)
		files[i]['PH_2'][0] = file_init['PH'][0] + PH_interval*(i+1)
		files[i]['QVAPOR'][0] = file_init['QVAPOR'][0] + QVAPOR_interval*(i+1)
		files[i]['Times'][0] = t
		files[i].close()

	file_init.close()

	return restarts






################################
#                              #
# Create an initial trajectory #
# that is a physical path      #
# according to the wrf model   #
#                              #
################################

def physical_profile(Nt, path_forward, path_input, storm_name, opt_interval, starts, run_storm):

	if run_storm == 'true':	
		shutil.copy2(path_input + 'namelist.input.forward_long', path_forward + 'namelist.input')
		os.chdir(path_forward)
		os.system('mpirun -np 28 ./wrf.exe')

		for filename in glob.glob(os.path.join(path_input, 'wrfout_d01*')):
			shutil.move(filename,  path_input + 'wrfout_' + storm_name + '.nc')

	else:
		shutil.copy2(path_input + 'namelist.input.forward_long_abbr', path_forward + 'namelist.input')
		os.chdir(path_forward)
		os.system('mpirun -np 28 ./wrf.exe')

	restarts = []
	
	for i in range(Nt+1):

		if i == 0:
			d1 = starts[2]/10
			d2 = starts[2]%10

			h1 = starts[3]/10
			h2 = starts[3]%10			
		else:
			h2 += opt_interval
			if h2 >= 10:
				h1 += h2/10
				h2 = h2%10
			if (h1 == 2 and h2 == 4):
				h1 = 0
				h2 = 0
				d2 += 1
				if d2%10 == 0:
					d1 += 1 
					d2 = 0
			
			cur_day = str(d1) + str(d2)
			cur_hour = str(h1) + str(h2)
			if starts[1] < 10:
				start_month = '0' + str(starts[1])
			else:
				start_month = str(starts[1])
				
			string = 'wrfrst_d01_' + str(starts[0]) + '-' + start_month + '-' + cur_day + '_' + cur_hour + ':00:00'
			restarts.append(string)
			shutil.move(path_input + string, path_input + string + '.nc')
			shutil.copy2(path_input + string + '.nc', path_input + 'ORIG_' + string + '.nc')

	return restarts








################################
#                              #
# Create an initial trajectory #
# that is a physical path      #
# according to the wrf model   # 
# plus a small perturbation    #
#                              #
################################

def perturb_profile(Nt, path_forward, path_input, storm_name, opt_interval, starts, run_storm):

	if run_storm == 'true':	
		shutil.copy2(path_input + 'namelist.input.forward_long', path_forward + 'namelist.input')
		os.chdir(path_forward)
		os.system('mpirun -np 28 ./wrf.exe')

		for filename in glob.glob(os.path.join(path_input, 'wrfout_d01*')):
			shutil.move(filename,  path_input + 'wrfout_' + storm_name + '.nc')

	else:
		shutil.copy2(path_input + 'namelist.input.forward_long_abbr', path_forward + 'namelist.input')
		os.chdir(path_forward)
		os.system('mpirun -np 28 ./wrf.exe')

	restarts = []
	
	for i in range(Nt+1):

		if i == 0:
			d1 = starts[2]/10
			d2 = starts[2]%10

			h1 = starts[3]/10
			h2 = starts[3]%10			
		else:
			h2 += opt_interval
			if h2 >= 10:
				h1 += h2/10
				h2 = h2%10
			if (h1 == 2 and h2 == 4):
				h1 = 0
				h2 = 0
				d2 += 1
				if d2%10 == 0:
					d1 += 1 
					d2 = 0
			
			cur_day = str(d1) + str(d2)
			cur_hour = str(h1) + str(h2)
			if starts[1] < 10:
				start_month = '0' + str(starts[1])
			else:
				start_month = str(starts[1])
				
			string = 'wrfrst_d01_' + str(starts[0]) + '-' + start_month + '-' + cur_day + '_' + cur_hour + ':00:00'
			restarts.append(string)
			shutil.move(path_input + string, path_input + string + '.nc')
			shutil.copy2(path_input + string + '.nc', path_input + 'ORIG_' + string + '.nc')

			dset = netCDF4.Dataset(path_input + string + '.nc', 'r+')
			P = dset['P'][:]
			for ii in range(np.shape(P)[0]):
				for jj in range(np.shape(P)[1]):
					for kk in range(np.shape(P)[2]):
						for ll in range(np.shape(P)[3]):
							P[ii,jj,kk,ll] += .05*(np.random.random()-.025)*P[ii,jj,kk,ll]
			dset['P'][:] = P

	return restarts



		



############################
#                          #
# Create namelists for     #
# forward and adjoint runs #
#                          #
############################

def get_namelists(Nt, path, ext, opt_interval):

	n1 = open(path + 'namelist.input.' +ext + '_0')
	n1_lines = n1.readlines()
	n1_lines_orig = copy.copy(n1_lines)
	i = 0
	stop = 0
	stm = 0
	run_hours = int(opt_interval)
	
	while (stop != 1):
		
		line = n1_lines[i].split()
		#if (line[0] == 'run_hours'):
		#	line2 = line[2].split(',')
		#	run_hours = int(line2[0])
		if (line[0] == 'start_hour'):
			line2 = line[2].split(',')
			sth = int(line2[0])
		if (line[0] == 'start_day'):
			line2 = line[2].split(',')
			std = int(line2[0])
		if (line[0] == 'end_day'):
			line2 = line[2].split(',')
			end = int(line2[0])
			stop = 1
		else:
			i += 1
			
	for i in range(1,Nt):
		
		sth += run_hours
		if (sth < (24 - run_hours)):
			enh = sth + run_hours
		elif (sth == 24):
			sth = 0
			std += 1
			enh = run_hours
			#end += 1
		else:
			enh = (sth + run_hours)%24
			end += 1
			
		f = open(path + 'namelist.input.' + ext + '_' + str(i),'w')
		
		for j in range(len(n1_lines)):
			
			line = n1_lines[j].split()
			if (len(line) >0 and line[0] == 'restart'):
				f.write('restart = .true.,' + '\n')			
			elif (len(line) >0 and line[0] == 'start_minute'):
				f.write('start_minute = ' +str(stm) +',' + '\n')
			elif (len(line) >0 and line[0] == 'start_hour'):
				f.write('start_hour = ' +str(sth) +',' + '\n')
                        elif (len(line) >0 and line[0] == 'start_day'):
				f.write('start_day = ' +str(std) +',' + '\n')
			elif (len(line) >0 and line[0] == 'end_hour'):
				f.write('end_hour = ' +str(enh) +',' + '\n')
			elif (len(line) >0 and line[0] == 'end_day'):
				f.write('end_day = ' +str(end) +',' + '\n')
			elif (len(line) >0 and line[0] == 'history_interval'):
				f.write('history_interval = ' +str(int(60*opt_interval)) +',' + '\n')
			else:
				f.write(n1_lines[j])


####################################
#                                  #
# Initialize the final sensitivity # 
# files for adjoint integration    #
#                                  #
####################################

def final_sens_init(Nt, path, opt_interval, starts):

	files = []
	P_big = []

	for i in range(Nt):

		files.append(netCDF4.Dataset(path + 'final_sens_d01_' + str(i) + '.nc','r+'))
		t = files[i]['Times'][:]

		if i == 0:
			d1 = starts[2]/10
			d2 = starts[2]%10

			h1 = (starts[3] + opt_interval)/10
			h2 = (starts[3] + opt_interval)%10
		else:
			h2 += opt_interval
			if h2 >= 10:
				h1 += h2/10
				h2 = h2%10
			if (h1 == 2 and h2 == 4):
				h1 = 0
				h2 = 0
				d2 += 1
				if d2%10 == 0:
					d1 += 1 
					d2 = 0

		t[0][8] = str(d1)
		t[0][9] = str(d2)
				
		t[0][11] = str(h1)
		t[0][12] = str(h2)

		files[i]['Times'][:] = t
		files[i].close()





##########################################
#                                        #
# Initialize input files, restart files  # 
# boundary conditions, final sensitivity #
# files, and namelists for optimization  #
#                                        #
##########################################

def initialize(dirs, params, storm_name, opt_interval, starts, phys_test, run_storm):

	path_input = dirs[1]
	path_forward = dirs[2]
	path_adj = dirs[3]
	path_output = dirs[4]
	path_met = dirs[5]
	Nt = params[0]
	alpha_min = params[3]
	R = params[4]

	for filename in glob.glob(os.path.join(path_met, 'met_em*')):
		shutil.copy2(filename, path_adj)
		shutil.copy2(filename, path_forward)

	shutil.copy2(path_input + 'namelist.input.create_IC', path_adj + 'namelist.input')
	os.chdir(path_adj)

	utils.log_write(path_output, 'Calling real.exe for initial conditions', 'true')

	os.system(path_adj + 'real.exe')                                         # create the input file
	shutil.copy2(path_adj + 'wrfinput_d01', path_input + 'wrfinput_d01_IC.nc')   # move initial condition to input directory

	utils.log_write(path_output, 'Calling real.exe for boundary conditions')

	shutil.copy2(path_input + 'namelist.input.forward_long', path_forward + 'namelist.input') 
	os.chdir(path_forward)
	os.system(path_forward + 'real.exe')                                     # create long run boundary conditions
	shutil.copy2(path_forward + 'wrfbdy_d01', path_adj + 'wrfbdy_d01')
	shutil.copy2(path_forward + 'wrfinput_d01', path_input + 'wrfinput_d01')

	shutil.copy2(path_input + 'namelist.input.forward', path_input + 'namelist.input.forward_0')
	shutil.copy2(path_input + 'namelist.input.adj', path_input + 'namelist.input.adj_0')

	utils.log_write(path_output, 'Creating namelists')

	get_namelists(Nt, path_input, 'forward', opt_interval)
	get_namelists(Nt, path_input, 'adj', opt_interval)

	utils.log_write(path_output , 'Creating initial profile')

	if phys_test == 'true':
		rsts = physical_profile(Nt, path_forward, path_input, storm_name, opt_interval, starts, run_storm)
	elif phys_test == 'pert':
		rsts = perturb_profile(Nt, path_forward, path_input, storm_name, opt_interval, starts, run_storm)
	else:
		rsts = initial_profile(Nt, path_forward, path_input, path_output, storm_name, opt_interval, starts, run_storm)

	utils.log_write(path_output, 'Initializing sensitivity files')

	shutil.copy2(path_input + 'namelist.input.forward_long', path_adj + 'namelist.input') 
	os.chdir(path_adj)
	os.system(path_adj + 'real.exe')  
	shutil.copy2(path_adj + 'wrfinput_d01', path_input + 'wrfinput_final_sens.nc')

	nc_copy.copy('wrfinput_final_sens.nc','f.nc', path_input)
	shutil.copy2(path_input + 'f.nc', path_input + 'final_sens_d01')
	for i in range(Nt+1):
		shutil.copy2(path_input + 'final_sens_d01', path_input + 'final_sens_d01_' + str(i) + '.nc')

	final_sens_init(Nt, path_input, opt_interval, starts)

	for filename in glob.glob(os.path.join(path_forward, 'met_em*')):
		os.remove(filename)
	for filename in glob.glob(os.path.join(path_adj, 'met_em*')):
		os.remove(filename)

	return rsts

	


