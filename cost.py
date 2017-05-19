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
import utils

#
# This module contains the routine for getting the cost funtion:
#
# 1) get_g: calculate the cost function from the 
#           perturbations and the final state
#

#########################
#                       #
# Get the cost function #                       
#                       #
#########################

#
# Nt                   Number of optimization timesteps
# g                    Total cost
# running_cost         Cost associated with perturbations
# final_cost           Cost associated with final state
# R                    Weight of final cost relative to running cost
# U_list               List of all perturbation variables for all timesteps
# U_P_full             List of P perturbations for all timesteps
# U_U_full             List of U perturbations for all timesteps
# U_V_full             List of V perturbations for all timesteps
# U_W_full             List of W perturbations for all timesteps
# U_MU_full            List of MU perturbations for all timesteps
# U_PH_full            List of PH perturbations for all timesteps
# U_QVAPOR_full        List of QVAPOR perturbations for all timesteps
# p_cont               Contribution of P to running cost
# u_cont               Contribution of U to running cost
# v_cont               Contribution of V to running cost
# w_cont               Contribution of W to running cost
# mu_cont              Contribution of MU to running cost
# ph_cont              Contribution of PH to running cost
# qvapor_cont          Contribution of QVAPOR to running cost
# p_final_cont         Contribution of P to running cost
# u_final_cont         Contribution of U to running cost
# v_final_cont         Contribution of V to running cost
# w_final_cont         Contribution of W to running cost
# mu_final_cont        Contribution of MU to running cost
# ph_final_cont        Contribution of PH to running cost
# qvapor_final_cont    Contribution of QVAPOR to running cost
# P_scale              Variance for P in cost function
# U_scale              Variance for U in cost function
# V_scale              Variance for V in cost function
# W_scale              Variance for W in cost function
# PH_scale             Variance for PH in cost function
# MU_scale             Variance for MU in cost function
# QVAPOR_scale         Variance for QVAPOR in cost function
#

def get_g(Nt, R, scales, U_list, dirs, restarts, storm_name):

	path_input = dirs[1]
	path_output = dirs[4]

	final_input = netCDF4.Dataset(path_input + restarts[Nt-1] + '.nc','r+')
	#final_output = netCDF4.Dataset(path_output + 'wrfout_' + str(Nt-1) + '.nc','r+')
	file_end = netCDF4.Dataset(path_input + 'wrfout_' + storm_name + '.nc','r+')

	P_end = file_end['P'][-1,:]
	U_end = file_end['U'][-1,:]
	V_end = file_end['V'][-1,:]
	W_end = file_end['W'][-1,:]
	MU_end = file_end['MU'][-1,:]
	PH_end = file_end['PH'][-1,:]
	QVAPOR_end = file_end['QVAPOR'][-1,:]

	utils.log_write(path_output, 'calling get_g')

	P_scale = scales[0]
	U_scale = scales[1]
	V_scale = scales[2]
	W_scale = scales[3]
	PH_scale = scales[4]
	MU_scale = scales[5]
	QVAPOR_scale = scales[6]

	U_P_full = U_list[0]
	U_U_full = U_list[1]
	U_V_full = U_list[2]
	U_W_full = U_list[3]
	U_PH_full = U_list[4]
	U_MU_full = U_list[5]
	U_QVAPOR_full = U_list[6]

	p_cont = (np.linalg.norm(U_P_full))**2
	u_cont = (np.linalg.norm(U_U_full))**2
	v_cont = (np.linalg.norm(U_V_full))**2
	w_cont = (np.linalg.norm(U_W_full))**2
	ph_cont = (np.linalg.norm(U_PH_full))**2
	mu_cont = (np.linalg.norm(U_MU_full))**2
	qvapor_cont = (np.linalg.norm(U_QVAPOR_full))**2

	p_final_cont = np.linalg.norm(final_input['P'][-1] - file_end['P'][70])**2 * 1./(P_scale**2)
	u_final_cont = np.linalg.norm(final_input['U_2'][-1] - file_end['U'][70])**2 * 1./(U_scale**2)
	v_final_cont = np.linalg.norm(final_input['V_2'][-1] - file_end['V'][70])**2 * 1./(V_scale**2)
	w_final_cont = np.linalg.norm(final_input['W_2'][-1] - file_end['W'][70])**2 * 1./(W_scale**2)
	ph_final_cont = np.linalg.norm(final_input['PH_2'][-1] - file_end['PH'][70])**2 * 1./(PH_scale**2)
	mu_final_cont = np.linalg.norm(final_input['MU_2'][-1] - file_end['MU'][70])**2 * 1./(MU_scale**2)
	qvapor_final_cont = np.linalg.norm(final_input['QVAPOR'][-1] - file_end['QVAPOR'][70])**2 * 1./(QVAPOR_scale**2)

	#p_final_cont = np.linalg.norm(final_output['P'][-1,:] - file_end['P'][70,:])**2 * 1./P_scale
	#u_final_cont = np.linalg.norm(final_output['U'][-1,:] - file_end['U'][70,:])**2 * 1./U_scale
	#v_final_cont = np.linalg.norm(final_output['V'][-1,:] - file_end['V'][70,:])**2 * 1./V_scale
	#w_final_cont = np.linalg.norm(final_output['W'][-1,:] - file_end['W'][70,:])**2 * 1./W_scale
	#ph_final_cont = np.linalg.norm(final_output['PH'][-1,:] - file_end['PH'][70,:])**2 * 1./PH_scale
	#mu_final_cont = np.linalg.norm(final_output['MU'][-1,:] - file_end['MU'][70,:])**2 * 1./MU_scale
	#qvapor_final_cont = np.linalg.norm(final_output['QVAPOR'][-1,:] - file_end['QVAPOR'][70,:])**2 * 1./QVAPOR_scale

#	p_final_cont = P_scale * np.linalg.norm(final_input['P'][:] - final_output['P'][-1,:])
#	u_final_cont = U_scale * np.linalg.norm(final_input['U_2'][:] - final_output['U'][-1,:])
#	v_final_cont = V_scale * np.linalg.norm(final_input['V_2'][:] - final_output['V'][-1,:])
#	w_final_cont = W_scale * np.linalg.norm(final_input['W_2'][:] - final_output['W'][-1,:])
#	ph_final_cont = PH_scale * np.linalg.norm(final_input['PH_2'][:] - final_output['PH'][-1,:])
#	mu_final_cont = MU_scale * np.linalg.norm(final_input['MU_2'][:] - final_output['MU'][-1,:])
#	qvapor_final_cont = QVAPOR_scale * np.linalg.norm(final_input['QVAPOR'][:] - final_output['QVAPOR'][-1,:])
	
	final_cost = p_final_cont + u_final_cont + v_final_cont + w_final_cont + ph_final_cont + mu_final_cont + qvapor_final_cont
	running_cost = p_cont + u_cont + v_cont + w_cont + ph_cont + mu_cont + qvapor_cont
	final_cost = final_cost/(R**2.)
	g = running_cost + final_cost

	utils.log_write(path_output, 'running cost is ' + str(running_cost) + ', final cost is ' + str(final_cost))
	utils.log_write(path_output, 'p cont is ' +str(p_cont) + ', u cont is ' +str(u_cont) + ', v cont is ' +str(v_cont) + ', w cont is ' +str(w_cont) + ', mu cont is ' +str(mu_cont) + ', ph cont is ' +str(ph_cont)+ ', qvapor cont is ' +str(qvapor_cont))
	utils.log_write(path_output, 'p final cont is ' +str(p_final_cont/R**2) + ', u final cont is ' +str(u_final_cont/R**2) + ', v final cont is ' +str(v_final_cont/R**2) + ', w final cont is ' +str(w_final_cont/R**2) + ', mu final cont is ' +str(mu_final_cont/R**2) + ', ph final cont is ' +str(ph_final_cont/R**2) + ', qvapor final cont is ' +str(qvapor_final_cont/R**2))
	utils.log_write(path_output, 'cost is ' + str(g))

	final_input.close()
	#final_output.close()
	file_end.close()
	
	return g


