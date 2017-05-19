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
# This module contains the routine creating a netCDF:
# file containing the final trajectory resulting
# from the minimization:
#
# 1) post_process
#

#######################################
#                                     #
# Puts final trajectory into one file #
#                                     #
#######################################

def post_process(params, dirs, rsts, opt_interval):

	path_input = dirs[1]
	path_output = dirs[4]

	utils.log_write(path_output, 'Calling post_process')

	Nt = params[0]
	dx = params[1]
	dy = params[2]
	files = []
	u_list = []
	v_list = []
	p_list = []
	q_list = []
	t = np.zeros(Nt+1)

	for i in range(Nt+1):

		if i == 0:
			files.append(netCDF4.Dataset(path_input + 'wrfinput_d01_IC.nc'))
			u_list.append(files[i]['U'][0, :, :, :-1])
			v_list.append(files[i]['V'][0, :, :-1, :])
			p_list.append(files[i]['P'][0])

		else:
			files.append(netCDF4.Dataset(path_input + 'OLD_' + rsts[i-1] + '.nc'))
			u_list.append(files[i]['U_1'][0, :, :, :-1])
			v_list.append(files[i]['V_1'][0, :, :-1, :])
			p_list.append(files[i]['P'][0])

		log_full = open(path_output + 'LOGFILE.txt','a')	
		log_full.write('timestep is ' + str(i) + ' ')
		log_full.write('max U is ' + str(np.max(u_list[i])))
		log_full.write('\n')
		log_full.close()

		t[i] = i*opt_interval

		q_cur = np.zeros(np.shape(u_list[i]))
		q_cur[:,:-1,:] = -1/dy*(u_list[i][:,1:,:] - u_list[i][:,:-1,:])
		q_cur[:,:,:-1] += 1/dx*(v_list[i][:,:,1:] - v_list[i][:,:,:-1])

		q_list.append(q_cur)

		if (i == 0):
			lat_old = files[i]['XLAT'][:]
			lon_old = files[i]['XLONG'][:]

		files[i].close()

	u_array = np.array(u_list)
	v_array = np.array(v_list)
	p_array = np.array(p_list)
	q_array = np.array(q_list)

	for i in range(len(u_array)):
		log_full = open(path_output + 'LOGFILE.txt','a')	
		log_full.write('max U is ' + str(np.max(u_array[i])))
		log_full.close()

	full_file = nc.netcdf_file(path_output + 'full_output.nc', 'w')
	full_file.createDimension('time', Nt+1)
	full_file.createDimension('pres', np.shape(u_array)[1])
	full_file.createDimension('lat', np.shape(u_array)[2])
	full_file.createDimension('lon', np.shape(u_array)[3])

	time = full_file.createVariable('time', 'float', ('time',))
	time[:] = t
	time.units = 'hours'

	lat = full_file.createVariable('lat', 'float', ('lat',))
	lat[:] = lat_old[0,:,0]
	lat.units = 'degrees_north'

	lon = full_file.createVariable('lon', 'float',  ('lon',))
	lon[:] = lon_old[0,0,:]
	lon.units = 'degrees_east'

	q = full_file.createVariable('q', 'float', ('time','pres','lat','lon'))
	q[:] = q_array
	q.units = 'm/s'

	p = full_file.createVariable('p', 'float', ('time','pres','lat','lon'))
	p[:] = p_array
	p.units = 'm/s'

	u = full_file.createVariable('u', 'float', ('time','pres','lat','lon'))
	u[:] = u_array
	u.units = 'm/s'

	v = full_file.createVariable('v', 'float', ('time','pres','lat','lon'))
	v[:] = v_array
	v.units = 'm/s'

