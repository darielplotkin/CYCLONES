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
import cost
import variables
import utils

#
# This module contains the routines for updating U variables
# at each interation of the minimization, as well as updating
# the final sensitivity files for adjoint integration.
# Routines are:
#
# 1) update_U: update the perturbation variables.
#
# 2) final_sens: update the sensitivity files.
#

##########################
#                        #
# Update the trajectory  #
# at each iteration step #                        
#                        #
##########################

#
# Nt                   Number of optimization timesteps
# alpha                Step size of the optimization
# alpha_min            Minimum step for the optimization
# R                    Weight associated with the final cost vs running cost (cost = running + final/R^2)
# p_P_old_list         List of conjugates from previous timestep for pressure
# p_U_old_list         List of conjugates from previous timestep for u
# p_V_old_list         List of conjugates from previous timestep for v
# p_W_old_list         List of conjugates from previous timestep for w
# p_MU_old_list        List of conjugates from previous timestep for mu
# grad_P_old_list      List of conjugates from previous timestep for pressure
# grad_U_old_list      List of conjugates from previous timestep for u
# grad_V_old_list      List of conjugates from previous timestep for v
# grad_W_old_list      List of conjugates from previous timestep for w
# grad_MU_old_list     List of conjugates from previous timestep for mu
# P_scale              Variance for P in cost function
# U_scale              Variance for U in cost function
# V_scale              Variance for V in cost function
# W_scale              Variance for W in cost function
# PH_scale             Variance for PH in cost function
# MU_scale             Variance for MU in cost function
# QVAPOR_scale         Variance for QVAPOR in cost function
# g_old                Cost function before update
# g_old                Cost function after update

def update_U(params, scales, restarts, storm_name, U_list, dirs, count, p_old_list = [], grad_old_list = [], IC_fixed = 'TRUE'):

############ get some parameters ##############
	Nt = params[0]
	alpha_min = params[3]
	R = params[4]
	path_input = dirs[1]
	path_forward = dirs[2]
	path_adj = dirs[3]
	path_output = dirs[4]
	P_scale = scales[0]
	U_scale = scales[1]
	V_scale = scales[2]
	W_scale = scales[3]
	PH_scale = scales[4]
	MU_scale = scales[5]
	QVAPOR_scale = scales[6]
	alpha = 1.

	utils.log_write(path_output, 'calling update_U, count is ' + str(count) )

############ if we're doing conjugate gradient, get the previous conjugates and gradients ############
	if len(p_old_list) > 0:

		p_P_old_list = p_old_list[0]
		p_U_old_list = p_old_list[1]
		p_V_old_list = p_old_list[2]
		p_W_old_list = p_old_list[3]
		grad_P_old_list = grad_old_list[0]
		grad_U_old_list = grad_old_list[1]
		grad_V_old_list = grad_old_list[2]
		grad_W_old_list = grad_old_list[3]

	else:

		p_P_old_list = []
		p_U_old_list = []
		p_V_old_list = []
		p_W_old_list = []
		grad_P_old_list = []
		grad_U_old_list = []
		grad_V_old_list = []
		grad_W_old_list = []

################ Copy the previous perturbations and states in case update leads to higher cost ##################
	U_list_OLD = copy.deepcopy(U_list)
	for i in range(Nt):
		shutil.copy2(path_input + restarts[i] + '.nc', path_input + 'OLD_' + restarts[i] + '.nc')

################ Get the old cost, initialize new cost to be larger than old cost ###################
	g_old = cost.get_g(Nt, R, scales, U_list, dirs, restarts, storm_name)
	g_new = 2*g_old

################## Create lists of restart files, output files, and file containing the final state ##############
	files_input = []
	files_output = []  # currently not used and probably not necessary
	file_end = netCDF4.Dataset(path_input + 'wrfout_' + storm_name + '.nc','r+')

	for i in range(Nt):

		files_input.append(netCDF4.Dataset(path_input + restarts[i] + '.nc','r+'))
		files_output.append(netCDF4.Dataset(path_output + 'wrfout_' + str(i) + '.nc','r+'))

################# Iterate over all timesteps #######################
	for i in range(Nt):


		j = Nt - i - 1     # j counts backwards from final time to initial time

		P_out = files_output[j]['P'][-1,:]  # getting the output of the j_th forward run, not currently used
		U_out = files_output[j]['U'][-1,:]
		V_out = files_output[j]['V'][-1,:]
		W_out = files_output[j]['W'][-1,:]
		MU_out = files_output[j]['MU'][-1,:]
		PH_out = files_output[j]['PH'][-1,:]
		QVAPOR_out = files_output[j]['QVAPOR'][-1,:]

		if i == 0:

################ if i = 0, then we initialize the adjoint variables as (target state - final state of current path) #########

			# final state of current path
			P_in = files_input[j]['P'][-1,:]       
			U_in = files_input[j]['U_2'][-1,:]
			V_in = files_input[j]['V_2'][-1,:]
			W_in = files_input[j]['W_2'][-1,:]
			MU_in = files_input[j]['MU_2'][-1,:]
			PH_in = files_input[j]['PH_2'][-1,:]
			QVAPOR_in = files_input[j]['QVAPOR'][-1,:]

			# target state
			P_end = file_end['P'][70,:]          
			U_end = file_end['U'][70,:]
			V_end = file_end['V'][70,:]
			W_end = file_end['W'][70,:]
			MU_end = file_end['MU'][70,:]
			PH_end = file_end['PH'][70,:]
			QVAPOR_end = file_end['QVAPOR'][70,:]

			# initialize the adjoint variables for final timestep
			P_adj = P_end - P_in               
			U_adj = U_end - U_in
			V_adj = V_end - V_in
			W_adj = W_end - W_in
			MU_adj = MU_end - MU_in
			PH_adj = PH_end - PH_in
			QVAPOR_adj = QVAPOR_end - QVAPOR_in

			max_adj = max(np.max(np.abs(U_adj)), np.max(np.abs(V_adj)))

			# create a list of all the adjoint variables
			adj_list = [P_adj, U_adj, V_adj, W_adj, MU_adj, PH_adj, QVAPOR_adj]

			# create final sensitivity file with the adjoint variables
			final_sens(j, path_input, adj_list, path_output)

############### running the adjoint model at the final time step #####################
			os.chdir(path_adj)
			# copy over the appropariate restart, namelist,  and final sensitivity files
			shutil.copy2(path_input + restarts[j-1] + '.nc', path_adj + restarts[j-1])
			shutil.copy2(path_input + 'final_sens_d01_' + str(j) + '.nc', path_adj + 'final_sens_d01')
			shutil.copy2(path_input + 'namelist.input.adj_' + str(j), path_adj + 'namelist.input')
			utils.log_write(path_output, 'calling adjoint model with timestep = ' +str(j))
			# call the model
			os.system('mpirun -np 28 ./wrf.exe')
			# remove old restart file
			os.remove(restarts[j-1])
			# move new restart file and log file to the output directory
			shutil.move(path_adj + 'rsl.out.0000', path_output + 'rsl.adj.' + str(count) + '.' + str(j))
			for filename in glob.glob(os.path.join(path_adj, 'gradient_wrfplus*')):
				shutil.move(filename, path_output + 'gradient_wrfplus_' + str(j) + '.nc')

			file_end.close()

		elif i < (Nt - 1):

			# read in adjoint variables from the gradient output from previous adjoint run
			file_gradient = netCDF4.Dataset(path_output + 'gradient_wrfplus_' + str(j+1) + '.nc','r+')

			P_adj = file_gradient['A_P'][:]
			U_adj = file_gradient['A_U'][:]
			V_adj = file_gradient['A_V'][:]
			W_adj = file_gradient['A_W'][:]
			MU_adj = file_gradient['A_MU'][:]
			PH_adj = file_gradient['A_PH'][:]
			QVAPOR_adj = file_gradient['A_QVAPOR'][:]

			file_gradient.close()

			# create a list of all the adjoint variables
			adj_list = [P_adj, U_adj, V_adj, W_adj, MU_adj, PH_adj, QVAPOR_adj]

			# create final sensitivity file with the adjoint variables
			final_sens(j, path_input, adj_list, path_output)

############### running the adjoint model at intermediate time steps #####################
			os.chdir(path_adj)
			# copy over the appropariate restart, namelist,  and final sensitivity files
			shutil.copy2(path_input + restarts[j-1] + '.nc', path_adj + restarts[j-1])
			shutil.copy2(path_input + 'final_sens_d01_' + str(j) + '.nc', path_adj + 'final_sens_d01')
			shutil.copy2(path_input + 'namelist.input.adj_' + str(j), path_adj + 'namelist.input')
			utils.log_write(path_output, 'calling adjoint model with timestep = ' +str(j))
			# call the model
			os.system('mpirun -np 28 ./wrf.exe')
			# remove old restart file
			os.remove(restarts[j-1])
			# move new restart file and log file to the output directory
			shutil.move(path_adj + 'rsl.out.0000', path_output + 'rsl.adj.' + str(count) + '.' + str(j))
			for filename in glob.glob(os.path.join(path_adj, 'gradient_wrfplus*')):
				shutil.move(filename, path_output + 'gradient_wrfplus_' + str(j) + '.nc')

		else:

 			# read in adjoint variables from the gradient output from previous adjoint run
			file_gradient = netCDF4.Dataset(path_output + 'gradient_wrfplus_' + str(j+1) + '.nc','r+')

			P_adj = file_gradient['A_P'][:]
			U_adj = file_gradient['A_U'][:]
			V_adj = file_gradient['A_V'][:]
			W_adj = file_gradient['A_W'][:]
			MU_adj = file_gradient['A_MU'][:]
			PH_adj = file_gradient['A_PH'][:]
			QVAPOR_adj = file_gradient['A_QVAPOR'][:]

			file_gradient.close()

			# create a list of all the adjoint variables
			adj_list = [P_adj, U_adj, V_adj, W_adj, MU_adj, PH_adj, QVAPOR_adj]

			# create final sensitivity file with the adjoint variables
			final_sens(j, path_input, adj_list, path_output)

############### running the adjoint model at the first time step #####################
			os.chdir(path_adj)
			# copy over the appropariate restart, namelist,  and final sensitivity files
			shutil.copy2(path_input + 'final_sens_d01_' + str(j) + '.nc', path_adj + 'final_sens_d01')
			shutil.copy2(path_input + 'namelist.input.adj_' + str(j), path_adj + 'namelist.input')
			utils.log_write(path_output, 'calling adjoint model with timestep = ' +str(j))
			# call the model
			os.system('mpirun -np 28 ./wrf.exe')
			# remove old restart file
			shutil.move(path_adj + 'rsl.out.0000', path_output + 'rsl.adj.' + str(count) + '.' + str(j))
			# move new restart file and log file to the output directory
			for filename in glob.glob(os.path.join(path_adj, 'gradient_wrfplus*')):
				shutil.move(filename, path_output + 'gradient_wrfplus_' + str(j) + '.nc')

	# initialize t
	P_update = copy.deepcopy(U_list[0])
	U_update = copy.deepcopy(U_list[1])
	V_update = copy.deepcopy(U_list[2])
	W_update = copy.deepcopy(U_list[3])
	PH_update = copy.deepcopy(U_list[4])
	MU_update = copy.deepcopy(U_list[5])
	QVAPOR_update = copy.deepcopy(U_list[6])

	files_gradient = []

	if (IC_fixed == 'TRUE'):
		start = 1
	else:
		start = 0

	for i in range(start, Nt):

		files_gradient.append(netCDF4.Dataset(path_output + 'gradient_wrfplus_' + str(i) + '.nc','r+'))

#	if (count == 1):
	if (count < 10000):

		for i in range(start, Nt):
	
			A_P = files_gradient[i-start]['A_P'][:]#/max(1., np.linalg.norm(files_gradient[i-start]['A_P']))
			A_U = files_gradient[i-start]['A_U'][:]#/max(1., np.linalg.norm(files_gradient[i-start]['A_U']))
			A_V = files_gradient[i-start]['A_V'][:]#/max(1., np.linalg.norm(files_gradient[i-start]['A_V']))
			A_W = files_gradient[i-start]['A_W'][:]#/max(1., np.linalg.norm(files_gradient[i-start]['A_W']))
			A_MU = files_gradient[i-start]['A_MU'][:]#/max(1., np.linalg.norm(files_gradient[i-start]['A_MU']))
			A_PH = files_gradient[i-start]['A_PH'][:]#/max(1., np.linalg.norm(files_gradient[i-start]['A_PH']))
			A_QVAPOR = files_gradient[i-start]['A_QVAPOR'][:]#/max(1., np.linalg.norm(files_gradient[i-start]['A_QVAPOR']))

			P_update[i] += P_scale * A_P/(R**2)
			U_update[i] += U_scale * A_U/(R**2)
			V_update[i] += V_scale * A_V/(R**2)
			W_update[i] += W_scale * A_W/(R**2)
			PH_update[i] += PH_scale * A_PH/(R**2)
			MU_update[i] += MU_scale * A_MU/(R**2)

#			grad_U_k1 = grad_U_k1/max(1,np.linalg.norm(A_U))#max(1,np.linalg.norm(grad_U_k1))
#			grad_V_k1 = grad_V_k1/max(1,np.linalg.norm(A_V))#max(1,np.linalg.norm(grad_V_k1))
#			grad_W_k1 = grad_W_k1/max(1,np.linalg.norm(A_W))#max(1,np.linalg.norm(grad_W_k1))
#			grad_P_k1 = grad_P_k1/max(1,np.linalg.norm(A_P))#max(1,np.linalg.norm(grad_P_k1))
#			grad_PH_k1 = grad_PH_k1/max(1,np.linalg.norm(A_PH))#max(1,np.linalg.norm(grad_PH_k1))
#			grad_MU_k1 = grad_MU_k1/max(1,np.linalg.norm(A_MU))#max(1,np.linalg.norm(grad_MU_k1))
#			grad_QVAPOR_k1 = grad_QVAPOR_k1/max(1,np.linalg.norm(A_QVAPOR))#max(1,np.linalg.norm(grad_QVAPOR_k1))

			p_P_old_list.append(-(P_update[i]))
			p_U_old_list.append(-(U_update[i]))
			p_V_old_list.append(-(U_update[i]))
			p_W_old_list.append(-(U_update[i]))
			grad_P_old_list.append((U_update[i]))
			grad_U_old_list.append((U_update[i]))
			grad_V_old_list.append((U_update[i]))
			grad_W_old_list.append((U_update[i]))
			
			files_gradient[i-start].close()

#	elif(count > 1):

#		for i in range(len(range(start, Nt))):

#			A_P = files_sens[i]['A_P'][:]
#			A_U = files_sens[i]['A_U'][:]
#			A_V = files_sens[i]['A_V'][:]
#			A_W = files_sens[i]['A_W'][:]
#			A_MU = files_sens[i]['A_MU'][:]
#			A_PH = files_sens[i]['A_PH'][:]
#			A_QVAPOR = files_sens[i]['A_QVAPOR'][:]

#			grad_U_k1 = files_gradient[i]['A_U'][:]
#			grad_V_k1 = files_gradient[i]['A_V'][:]
#			grad_W_k1 = files_gradient[i]['A_W'][:]
#			grad_P_k1 = files_gradient[i]['A_P'][:]
#			grad_PH_k1 = files_gradient[i]['A_PH'][:]
#			grad_MU_k1 = files_gradient[i]['A_MU'][:]
#			grad_QVAPOR_k1 = files_gradient[i]['A_QVAPOR'][:]

#			grad_U_k1 = grad_U_k1/max(1,np.linalg.norm(A_U))#max(1,np.linalg.norm(grad_U_k1))
#			grad_V_k1 = grad_V_k1/max(1,np.linalg.norm(A_V))#max(1,np.linalg.norm(grad_V_k1))
#			grad_W_k1 = grad_W_k1/max(1,np.linalg.norm(A_W))#max(1,np.linalg.norm(grad_W_k1))
#			grad_P_k1 = grad_P_k1/max(1,np.linalg.norm(A_P))#max(1,np.linalg.norm(grad_P_k1))
#			grad_PH_k1 = grad_PH_k1/max(1,np.linalg.norm(A_PH))#max(1,np.linalg.norm(grad_PH_k1))
#			grad_MU_k1 = grad_MU_k1/max(1,np.linalg.norm(A_MU))#max(1,np.linalg.norm(grad_MU_k1))
#			grad_QVAPOR_k1 = grad_QVAPOR_k1/max(1,np.linalg.norm(A_QVAPOR))#max(1,np.linalg.norm(grad_QVAPOR_k1))

#			grad_P_new = (P_diff_k - grad_P_k1)/P_scale
#			grad_U_new = (U_diff_k - grad_U_k1)/U_scale
#			grad_V_new = (V_diff_k - grad_V_k1)/V_scale
#			grad_W_new = (W_diff_k - grad_W_k1)/W_scale

#			beta_P = (np.linalg.norm(grad_P_new)/np.linalg.norm(grad_P_old_list[i-start]))**2
#			beta_U = (np.linalg.norm(grad_U_new)/np.linalg.norm(grad_U_old_list[i-start]))**2
#			beta_V = (np.linalg.norm(grad_V_new)/np.linalg.norm(grad_V_old_list[i-start]))**2
#			beta_W = (np.linalg.norm(grad_W_new)/np.linalg.norm(grad_W_old_list[i-start]))**2

#			print('Betas are: B_p = ' +str(beta_P) + ', B_u = ' +str(beta_U)+ ', B_v = ' +str(beta_V)+ ', B_w = ' +str(beta_W))

#			p_P_new = -grad_P_new + beta_P * p_P_old_list[i-start]
#			p_U_new = -grad_U_new + beta_U * p_U_old_list[i-start]
#			p_V_new = -grad_V_new + beta_V * p_V_old_list[i-start]
#			p_W_new = -grad_W_new + beta_W * p_W_old_list[i-start]

#			files_input[i]['P'][:] += alpha*p_P_new
#			files_input[i]['U_1'][:] += alpha*p_U_new
#			files_input[i]['V_1'][:] += alpha*p_V_new
#			files_input[i]['U_2'][:] += alpha*p_U_new
#			files_input[i]['V_2'][:] += alpha*p_V_new
#			files_input[i]['W_1'][:] += alpha*p_W_new
#			files_input[i]['W_2'][:] += alpha*p_W_new
			#files_input[i]['PH_1'][:] = PH_1_k + alpha*(PH_diff_k - grad_PH_k1)/PH_scale
			#files_input[i]['MU_1'][:] = MU_1_k + alpha*(MU_diff_k - grad_MU_k1)/MU_scale
			#files_input[i]['PH_2'][:] = PH_1_k + alpha*(PH_diff_k - grad_PH_k1)/PH_scale
			#files_input[i]['MU_2'][:] = MU_1_k + alpha*(MU_diff_k - grad_MU_k1)/MU_scale
			#files_input[i]['QVAPOR'][:] = QVAPOR_1_k + alpha*(QVAPOR_diff_k - grad_QVAPOR_k1)/QVAPOR_scale

#			p_P_old_list[i-start] = p_P_new
#			p_U_old_list[i-start] = p_U_new
#			p_V_old_list[i-start] = p_V_new
#			p_W_old_list[i-start] = p_W_new
#			grad_P_old_list[i-start] = grad_P_new
#			grad_U_old_list[i-start] = grad_U_new
#			grad_V_old_list[i-start] = grad_V_new
#			grad_W_old_list[i-start] = grad_W_new
			
#			files_gradient[i-start].close()

	else:

		print('Counter is fucked up in update_inputs')

	#max_update = max(np.max(np.abs(U_update[Nt-1])), np.max(np.abs(V_update[Nt-1])))
	#utils.log_write(path_output, 'max_update is ' + str(max_update))
#
#	if (alpha*max_update > 5):
#
#		scale_factor = (alpha*max_update)/5.
#		utils.log_write(path_output, 'scale_factor is ' + str(scale_factor))
#		alpha = alpha/scale_factor
	adj_diff = 10
	while(g_new >= g_old and alpha > alpha_min or adj_diff > 5):

		U_list[0] -= alpha*np.array(P_update)
		U_list[1] -= alpha*np.array(U_update) 
		U_list[2] -= alpha*np.array(V_update) 
		U_list[3] -= alpha*np.array(W_update) 
		U_list[4] -= alpha*np.array(PH_update) 
		U_list[5] -= alpha*np.array(MU_update) 
		U_list[6] -= alpha*np.array(QVAPOR_update)

		variables.U_to_x(Nt, restarts, U_list, dirs, count)
		
		file_end = netCDF4.Dataset(path_input + 'wrfout_' + storm_name + '.nc','r+')
		last_rst = netCDF4.Dataset(path_input + restarts[i] + '.nc','r+')
		last_U = file_end['U'][70,:]
		last_V = file_end['V'][70,:]
		new_U = last_rst['U_2'][-1,:]
		new_V = last_rst['V_2'][-1,:]
		new_adj_U = np.max(np.abs(last_U - new_U))
		new_adj_V = np.max(np.abs(last_V - new_V))
		new_max_adj = max(new_adj_U, new_adj_U)
		adj_diff = new_max_adj - max_adj
		utils.log_write(path_output, 'adj_diff is ' + str(adj_diff))

		g_new = cost.get_g(Nt, R, scales, U_list, dirs, restarts, storm_name)
		alpha = alpha/2.

		utils.log_write(path_output, 'alpha is ' + str(alpha) + '       g_new is ' + str(g_new) + '       g_old is ' + str(g_old))

		if (g_new >= g_old or adj_diff > 5):
			for i in range(Nt):
				shutil.copy2(path_input + 'OLD_' + restarts[i] + '.nc', path_input + restarts[i] + '.nc')
			U_list = copy.deepcopy(U_list_OLD)

	if alpha < alpha_min:
		
		utils.log_write(path_output, 'cannot find a descent direction')
		g_new = g_old

	p_old_list = [p_P_old_list, p_U_old_list, p_V_old_list, p_W_old_list]
	grad_old_list = [grad_P_old_list, grad_U_old_list, grad_V_old_list, grad_W_old_list]

	return p_old_list, grad_old_list, U_list, g_new







############################
#                          #
# Update final sensitivies #
# at each iteration step   #                        
#                          #
############################

def final_sens(j, path_input, adj_list, path_output):

	utils.log_write(path_output, 'calling final_sens with j = ' + str(j))
	file_final = netCDF4.Dataset(path_input + 'final_sens_d01_' + str(j) + '.nc','r+')

	if (np.shape(adj_list[0])[0] != 1):
		for i in range(len(adj_list)):
			adj_list[i] = np.expand_dims(adj_list[i], axis = 0)

	file_final['A_P'][:] = adj_list[0]
	file_final['A_U'][:] = adj_list[1]
	file_final['A_V'][:] = adj_list[2]
	file_final['A_W'][:] = adj_list[3]
	file_final['A_MU'][:] = adj_list[4]
	file_final['A_PH'][:] = adj_list[5]
	file_final['A_QVAPOR'][:] = adj_list[6]

	file_final.close()



