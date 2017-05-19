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
import utils

#
# This module contains routines for transforming
# between U (perturbation) and x (physical) variables.
# Routines are:
#
# 1) x_to_U
#
# 2) U_to_x
#

##################################
#                                #
# Get the perturbation variables #
# from the physical variables    #                        
#                                #
##################################

def x_to_U(restarts, Nt, scales, dirs):

	path_input = dirs[1]
	path_output = dirs[4]

	utils.log_write(path_output, 'calling x_to_U')

	P_scale = scales[0]
	U_scale = scales[1]
	V_scale = scales[2]
	W_scale = scales[3]
	PH_scale = scales[4]
	MU_scale = scales[5]
	QVAPOR_scale = scales[6]
	
	files_output = []
	files_input = []

	U_P_full = [] 
	U_U_full = [] 
	U_V_full = [] 
	U_W_full = []
	U_PH_full = []
	U_MU_full = []
	U_QVAPOR_full = [] 

	for i in range(Nt):

		files_input.append(netCDF4.Dataset(path_input + restarts[i] + '.nc','r+'))
		files_output.append(netCDF4.Dataset(path_output + 'wrfout_' + str(i) + '.nc','r+'))

	for i in range(Nt):
	
		P_out = files_output[i]['P'][-1,:]
		U_out = files_output[i]['U'][-1,:]
		V_out = files_output[i]['V'][-1,:]
		W_out = files_output[i]['W'][-1,:]
		MU_out = files_output[i]['MU'][-1,:]
		PH_out = files_output[i]['PH'][-1,:]
		QVAPOR_out = files_output[i]['QVAPOR'][-1,:]

		P_in = files_input[i]['P'][:]
		U_in = files_input[i]['U_2'][:]
		V_in = files_input[i]['V_2'][:]
		W_in = files_input[i]['W_2'][:]
		MU_in = files_input[i]['MU_2'][:]
		PH_in = files_input[i]['PH_2'][:]
		QVAPOR_in = files_input[i]['QVAPOR'][:]

		if (np.shape(P_out)[0] != 1):
			P_out = np.expand_dims(P_out, axis = 0)
			U_out = np.expand_dims(U_out, axis = 0)
			V_out = np.expand_dims(V_out, axis = 0)
			W_out = np.expand_dims(W_out, axis = 0)
			PH_out = np.expand_dims(PH_out, axis = 0)
			MU_out = np.expand_dims(MU_out, axis = 0)
			QVAPOR_out = np.expand_dims(QVAPOR_out, axis = 0)
			#QVAPOR_out_tmp = np.zeros([1, np.shape(QVAPOR_out)[0], np.shape(QVAPOR_out)[1], np.shape(QVAPOR_out)[2]])
			#QVAPOR_out_tmp[0] = QVAPOR_out
			#QVAPOR_out = QVAPOR_out_tmp

		U_P = (P_in - P_out)/P_scale
		U_U = (U_in - U_out)/U_scale
		U_V = (V_in - V_out)/V_scale
		U_W = (W_in - W_out)/W_scale
		U_MU = (MU_in - MU_out)/MU_scale
		U_QVAPOR = (QVAPOR_in - QVAPOR_out)/QVAPOR_scale
		U_PH = (PH_in - PH_out)/PH_scale
	
		U_P_full.append(U_P)
		U_U_full.append(U_U)
		U_V_full.append(U_V)
		U_W_full.append(U_W)
		U_PH_full.append(U_PH)
		U_MU_full.append(U_MU)
		U_QVAPOR_full.append(U_QVAPOR)

	U_P_full = np.array(U_P_full)
	U_U_full = np.array(U_U_full)
	U_V_full = np.array(U_V_full)
	U_W_full = np.array(U_W_full)
	U_MU_full = np.array(U_MU_full)
	U_QVAPOR_full = np.array(U_QVAPOR_full)
	U_PH_full = np.array(U_PH_full)

	utils.log_write(path_output, 'maxes of U_U are ' +str(np.max(U_U_full[0]))+"," +str(np.max(U_U_full[1]))+"," +str(np.max(U_U_full[2]))+"," +str(np.max(U_U_full[3])) )

	U_list = [U_P_full, U_U_full, U_V_full, U_W_full, U_PH_full, U_MU_full, U_QVAPOR_full]

	return U_list




###################################
#                                 #
# Get the physical variables from #
# the perturbation variables      #                        
#                                 #
###################################


def U_to_x(Nt, restarts, U_list, dirs, count):

	path_input = dirs[1]
	path_forward = dirs[2]
	path_adj = dirs[3]
	path_output = dirs[4]
	path_met = dirs[5]
	utils.log_write(path_output, 'calling U_to_x')

	U_P_full = U_list[0]
	U_U_full = U_list[1]
	U_V_full = U_list[2]
	U_W_full = U_list[3]
	U_PH_full = U_list[4]
	U_MU_full = U_list[5]
	U_QVAPOR_full = U_list[6]

#	for i in range(Nt):
#
#		file_rst = netCDF4.Dataset(path_input + restarts[i] + '.nc','r+')
#		file_rst['P'][:] += U_P_full[i]
#		file_rst['U_1'][:] += U_U_full[i]
#		file_rst['V_1'][:] += U_V_full[i]
#		file_rst['W_1'][:] += U_W_full[i]
#		file_rst['U_2'][:] += U_U_full[i]
#		file_rst['V_2'][:] += U_V_full[i]
#		file_rst['W_2'][:] += U_W_full[i]
#		file_rst['PH_1'][:] += U_PH_full[i]
#		file_rst['MU_1'][:] += U_MU_full[i]
#		file_rst['PH_2'][:] += U_PH_full[i]
#		file_rst['MU_2'][:] += U_MU_full[i]
#		file_rst['QVAPOR'][:] += U_QVAPOR_full[i]
#		utils.log_write(path_output, 'i is '+ str(i) + ', max U_2 is ' + str(np.max(file_rst['U_2'][:])))

	for i in range(Nt):

		if i == 0:

			os.chdir(path_forward)
			shutil.copy2(path_input + 'namelist.input.forward_0', path_forward + 'namelist.input')
			os.system('mpirun -np 28 ./wrf.exe')
			shutil.move(path_forward + 'rsl.out.0000', path_output + 'rsl.' + str(count) + '.' + str(i))
			for filename in glob.glob(os.path.join(path_forward, 'wrfout*')):
				shutil.move(filename, path_output + 'wrfout_' + str(i) + '.nc')
			shutil.move(path_forward + restarts[i], path_input + restarts[i] + '.nc')

		else: 

			os.chdir(path_forward)
			shutil.copy2(path_input + restarts[i-1] + '.nc', path_forward + restarts[i-1])
			shutil.copy2(path_input + 'namelist.input.forward_' + str(i), path_forward + 'namelist.input')
			os.system('mpirun -np 28 ./wrf.exe')

			os.remove(restarts[i-1])
			shutil.move(path_forward + 'rsl.out.0000', path_output + 'rsl.' + str(count) + '.' + str(i))
			for filename in glob.glob(os.path.join(path_forward, 'wrfout*')):
				shutil.move(filename, path_output + 'wrfout_' + str(i) + '.nc')
			shutil.move(path_forward + restarts[i], path_input + restarts[i] + '.nc')

		file_rst = netCDF4.Dataset(path_input + restarts[i] + '.nc','r+')
		file_rst['P'][:] += U_P_full[i]
		file_rst['U_1'][:] += U_U_full[i]
		file_rst['V_1'][:] += U_V_full[i]
		file_rst['W_1'][:] += U_W_full[i]
		file_rst['U_2'][:] += U_U_full[i]
		file_rst['V_2'][:] += U_V_full[i]
		file_rst['W_2'][:] += U_W_full[i]
		file_rst['PH_1'][:] += U_PH_full[i]
		file_rst['MU_1'][:] += U_MU_full[i]
		file_rst['PH_2'][:] += U_PH_full[i]
		file_rst['MU_2'][:] += U_MU_full[i]
		file_rst['QVAPOR'][:] += U_QVAPOR_full[i]
		utils.log_write(path_output, 'Now i is '+ str(i) + ', max U_2 is ' + str(np.max(file_rst['U_2'][:])))

		file_rst.close()

