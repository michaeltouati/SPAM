#######################################################################
##                                                                   ##
##  Stopping Power of Protons and Alpha particles in Ambient Matter  ##
##                              (SPAM)                               ##
##                                                                   ##
## Copyright © 2019 Michaël J TOUATI                                 ##
##                                                                   ##
## This file is part of SPAM.                                        ##
##                                                                   ##
## SPAM is free software: you can redistribute it and/or modify      ##
## it under the terms of the GNU General Public License as published ##
## by the Free Software Foundation, either version 3 of the License, ##
## or (at your option) any later version.                            ##
##                                                                   ##
## SPAM is distributed in the hope that it will be useful,           ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of    ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     ##
## GNU General Public License for more details.                      ##
##                                                                   ##
## You should have received a copy of the GNU General Public License ##
## along with ESVM. If not, see <https://www.gnu.org/licenses/>.     ##
##                                                                   ##
#######################################################################
import tkinter
from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image
import os
import os.path
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
import math
import physical_constants
import maths_library
import physics_library

# import array of elements
s  = physical_constants.symbol
# proton mass :
mp = physical_constants.proton_mass
# proton mass :
mp = physical_constants.proton_mass
# alpha particle mass
mu = physical_constants.atomic_mass_unit
Aa = physical_constants.standard_atomic_weight_He
ma = Aa * mu 

def close(*args):
	try :
	    sys.exit()
	except ValueError:
		print('ERROR IN FUNCTION : close')
		pass

def generate_material_list():
# input  : nothing
# output : table = list of all implemented materials
	try :
		N     = len(s)
		table = []
		for n in range(0,N) :
			mat = eval("physical_constants.name_"+s[n])
			if ( len(mat) > 10 ) : # no elements with more than 10 different phases
				Z = physics_library.get_atomic_number(mat)
				k = int(Z) - 1
				name = '['+s[k]+']'+str(Z)+' : '+ mat
				table.append(name)
			else :
				for m in range(0,len(mat)):
					mat2 = mat[m]
					Z    = physics_library.get_atomic_number(mat2)
					k    = int(Z) - 1
					name = '['+s[k]+']'+str(Z)+' : '+ mat2
					table.append(name)
		table.append('GafChromic HD-V2 (solid at room temperature)')
		table.append('Polyester (solid at room temperature)')
		table.append('Plastic scintillator (vinyltoluene-based solid at room temperature)')
		table.append('Water (liquid at room temperature)')
		table.append('Air (dry gas at room temperature)')
		table.append('Glass SiO2 (solid at room temperature)')
		return table;
	except ValueError:
		print('ERROR IN FUNCTION : generate_material_list')
		pass

def get_material(target):
# input  : target_var = material chosen by the user (string of characters)
# output : material   = material name (string of characters)
	try :
		N             = len(s)
		material = ''
		for n in range(0,N) :
			mat = eval("physical_constants.name_"+s[n])
			if ( len(mat) > 10 ) : # no elements with more than 10 different phases
				Z = physics_library.get_atomic_number(mat)
				k = int(Z) - 1
				name = '['+s[k]+']'+str(Z)+' : '+ mat
				if (target == name):
					material = mat
			elif ( len(mat) <= 9 ) :
				for m in range(0,len(mat)):
					mat2 = mat[m]
					Z    = physics_library.get_atomic_number(mat2)
					k    = int(Z) - 1
					name = '['+s[k]+']'+str(Z)+' : '+ mat2
					if (target == name):
						material = mat2
		if (material == '') : # For compounds
			if (target == 'GafChromic HD-V2 (solid at room temperature)'):
				material = target
			elif (target == 'Polyester (solid at room temperature)'):
				material = target
			elif (target =='Plastic scintillator (vinyltoluene-based solid at room temperature)'):
				material = target
			elif (target == 'Water (liquid at room temperature)'):
				material = target
			elif (target == 'Air (dry gas at room temperature)') :
				material = target
			elif (target == 'Glass SiO2 (solid at room temperature)') :
				material = target
		return material;
	except ValueError:
		print('ERROR IN FUNCTION : get_material')
		pass

def calculate_left(*args):
# left button action
    try:
    	projectile_value = str(projectile_entry.get())
    	target_value     = str(target_entry.get())
    	compute_plot_and_save_stopping_power(projectile_value,target_value)
    except ValueError:
    	print('ERROR IN FUNCTION : calculate_left')
    	pass

def calculate_right(*args):
# right button action
	try:
		projectile_value = str(projectile_entry.get())
		target_value     = str(target_entry.get())
		Eps_ini_value    = float(Eps_ini_entry.get())
		compute_plot_and_save_Bragg_peak(projectile_value,target_value,Eps_ini_value)
	except ValueError:
		print('ERROR IN FUNCTION : calculate_right')
		pass

def plot_stopping_power(projectile_var,target_var,Eps_var,S_var):
# input  : projectile_var = 'Proton' or 'Alpha particle'
#          target_var     = name of the target material (ex. : "Boron (solid at room temperature)")
#          Eps_var        = array containing all projectile kinetic energy in (MeV)
#          S_var          = array containing the corresponding value of the stopping power in (MeV/cm)
# output : nothing - it just plots the stopping power versus the 
#                    energy together with the corresponding NIST data (if it exists) 
#                and save the corresponding figure dE_ds.png
	try :
		# check if NIST benchmark exists
		material  = get_material(target_var)
		Z         = physics_library.get_atomic_number(material)
		filename_e  = ''
		filename_n  = ''
		if (int(Z) == 6) :
			if (physical_constants.name_C[0] == material):
				if ( str(projectile_var) == 'Proton' ) :
					filename_e = '/share/apps/spam-1.0/benchmark/pstar/electronic/C_amorpheous.dat'
					filename_n = '/share/apps/spam-1.0/benchmark/pstar/nuclear/C_amorpheous.dat'
				elif ( str(projectile_var) == 'Alpha particle' ) :
					filename_e = '/share/apps/spam-1.0/benchmark/astar/electronic/C_amorpheous.dat'
					filename_n = '/share/apps/spam-1.0/benchmark/astar/nuclear/C_amorpheous.dat'
			elif (physical_constants.name_C[1] == material):
				if ( str(projectile_var) == 'Proton' ) :
					filename_e = '/share/apps/spam-1.0/benchmark/pstar/electronic/C_graphite.dat'
					filename_n = '/share/apps/spam-1.0/benchmark/pstar/nuclear/C_graphite.dat'
				elif ( str(projectile_var) == 'Alpha particle' ) :
					filename_e = '/share/apps/spam-1.0/benchmark/astar/electronic/C_graphite.dat'
					filename_n = '/share/apps/spam-1.0/benchmark/astar/nuclear/C_graphite.dat'
		elif (int(Z) == 50) :
			if (physical_constants.name_Sn[0] == material) :
				if ( str(projectile_var) == 'Proton' ) :
					filename_e = str('/share/apps/spam-1.0/benchmark/pstar/electronic/Sn.dat')
					filename_n = str('/share/apps/spam-1.0/benchmark/pstar/nuclear/Sn.dat')
				elif ( str(projectile_var) == 'Alpha particle' ) :
					filename_e = str('/share/apps/spam-1.0/benchmark/astar/electronic/Sn.dat')
					filename_n = str('/share/apps/spam-1.0/benchmark/astar/nuclear/Sn.dat')
		elif isinstance(Z, int) :
			if ( str(projectile_var) == 'Proton' ) :
				filename_e  = str('/share/apps/spam-1.0/benchmark/pstar/electronic/'+s[int(Z)-1]+'.dat')
				filename_n  = str('/share/apps/spam-1.0/benchmark/pstar/nuclear/'+s[int(Z)-1]+'.dat')
			elif ( str(projectile_var) == 'Alpha particle' ) :
				filename_e  = str('/share/apps/spam-1.0/benchmark/astar/electronic/'+s[int(Z)-1]+'.dat')
				filename_n  = str('/share/apps/spam-1.0/benchmark/astar/nuclear/'+s[int(Z)-1]+'.dat')
		elif (material == 'Plastic scintillator (vinyltoluene-based solid at room temperature)') :
			if ( str(projectile_var) == 'Proton' ) :
				filename_e = '/share/apps/spam-1.0/benchmark/pstar/electronic/Plastic_scintillator.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/pstar/nuclear/Plastic_scintillator.dat'
			else :
				filename_e = '/share/apps/spam-1.0/benchmark/astar/electronic/Plastic_scintillator.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/astar/nuclear/Plastic_scintillator.dat'
		elif (material == 'Water (liquid at room temperature)'):
			if ( str(projectile_var) == 'Proton' ) :
				filename_e = '/share/apps/spam-1.0/benchmark/pstar/electronic/Water.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/pstar/nuclear/Water.dat'
			else :
				filename_e = '/share/apps/spam-1.0/benchmark/astar/electronic/Water.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/astar/nuclear/Water.dat'
		elif (material == 'Air (dry gas at room temperature)') :
			if ( str(projectile_var) == 'Proton' ) :
				filename_e = '/share/apps/spam-1.0/benchmark/pstar/electronic/Air.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/pstar/nuclear/Air.dat'
			else :
				filename_e = '/share/apps/spam-1.0/benchmark/astar/electronic/Air.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/astar/nuclear/Air.dat'
		elif (material == 'Glass SiO2 (solid at room temperature)') :
			if ( str(projectile_var) == 'Proton' ) :
				filename_e = '/share/apps/spam-1.0/benchmark/pstar/electronic/Glass.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/pstar/nuclear/Glass.dat'
			else :
				filename_e = '/share/apps/spam-1.0/benchmark/astar/electronic/Glass.dat'
				filename_n = '/share/apps/spam-1.0/benchmark/astar/nuclear/Glass.dat'
		else :
			filename_e  = ''
			filename_n  = ''
		# create image dE_ds.png
		ymin      = min(S_var)
		n_ymin    = int(math.log(ymin)/math.log(10.)) - 1
		ymax      = max(S_var)
		n_ymax    = 2 + int(math.log(ymax)/math.log(10.))
		ymin_plot = float(10.**n_ymin)
		ymax_plot = float(10.**n_ymax)
		xmin_plot = 1.e-3
		if ( str(projectile_var) == 'Proton' ) :
			xmax_plot = 1.e4
		elif ( str(projectile_var) == 'Alpha particle' ) :
			xmax_plot = 1.e5
		font = {'style':  'normal',
		'color':  'black',
		'weight': 'normal',
		'size': 16,
				}
		fig=plt.figure()
		plt.loglog(Eps_var,S_var ,'red',linestyle='-' ,linewidth=2,label='SPAM')
		plt.grid(True,which="both",axis="both",color='gray', linestyle='-', linewidth=1)
		plt.xlim([xmin_plot,xmax_plot])
		plt.xticks(fontsize=16)
		plt.xlabel('E (MeV)', fontdict=font)
		plt.ylabel('S (MeV/cm)', fontdict=font)
		plt.yticks(fontsize=16)
		plt.ylim([ymin_plot,ymax_plot])
		ttl = str(projectile_var)+' Stopping Power in \n'+str(target_var)
		plt.title(ttl, fontdict=font)
		Eps_NIST = []
		S_NIST   = []
		if os.path.exists(filename_e) == True :
			rho  = physics_library.get_density(material)
			file = open(filename_e, 'r')
			for line in file:
				line      = line.strip()
				array     = line.split()
				Eps_NIST.append(float(array[0]))   
				S_NIST.append(float(array[1])*rho)
			file.close()
			file = open(filename_n, 'r')
			i = 0
			for line in file:
				line        = line.strip()
				array       = line.split()  
				S_NIST[i]   = S_NIST[i] + (float(array[1])*rho)
				i           = i + 1
			file.close()
			plt.loglog(Eps_NIST,S_NIST,'black',linestyle='--' ,linewidth=2,label='NIST')
			leg = plt.legend(loc='best', fontsize=16, fancybox=True)
			leg.get_frame().set_alpha(0.5)
		fig.savefig("dE_ds.png",bbox_inches='tight')
		print('* Stopping power plot is saved')
	except ValueError:
		print('ERROR IN FUNCTION : plot_stopping_power')
		pass

def save_stopping_power(N_gridpoints_var,Eps_var,S_var) :
# input  : N_gridpoints   = number of gridpoints
#          Eps_var        = array containing all projectile kinetic energy in (MeV)
#          S_var          = array containing the corresponding value of the stopping power in (MeV/cm)
# output : nothing - it just saves the text file dE_ds.txt containing the stopping power versus the 
#                    energy together
	try :
		# create dE_ds.txt file
		txtfile = open("dE_ds.txt","w+")
		txtfile.write('      E (MeV)     '+' '+'   S (MeV/cm)    '+" \n")
		for i in range(0,N_gridpoints_var):
			Eps_str = "{0:.15f}".format(Eps_var[i]) 
			S_str   = "{0:.15f}".format(S_var[i])
			txtfile.write(Eps_str+' '+S_str+" \n")
		txtfile.close()
		print('* Stopping power text file is saved')
		print('-----------------------------------')
	except ValueError:
		print('ERROR IN FUNCTION : save_stopping_power')
		pass

def compute_plot_and_save_stopping_power(projectile_var,target_var):
# input  : projectile_var = 'Proton' or 'Alpha particle'
#          target_var     = name of the target material (ex. : "Boron (solid at room temperature)")
# output : nothing - it just plots the projectile stopping power versus its kinetic energy
#                    save the corresponding figure dE_ds.png, 
#                    and generate the corresponding text file dE_ds.txt containing all values 
	try:
		# projectile
		if ( str(projectile_var) == 'Proton' ) :
			z = 1
			m = mp
		elif ( str(projectile_var) == 'Alpha particle' ) :
			z = 2
			m = ma
		# target
		material = get_material(target_var)
		# Stopping power versus Kinetic energy
		print('* Compute Stopping power curve')
		n_gridpoints = np.arange(-3.,5.,0.01)
		N_gridpoints = len(n_gridpoints)
		Eps          = np.zeros(N_gridpoints)
		S            = np.zeros(N_gridpoints)
		for i in range(0,N_gridpoints):
			Eps[i] = 10.**n_gridpoints[i]
			S[i]   = physics_library.Stopping_power(material,z,m,Eps[i])
		# plot
		plot_stopping_power(projectile_var,target_var,Eps,S)
		# save
		save_stopping_power(N_gridpoints,Eps,S)
		# display
		img = Image.open('dE_ds.png')
		img = img.resize((460, 355), Image.ANTIALIAS)
		img = ImageTk.PhotoImage(img)
		ttk.Label(mainframe, image = img).grid(column=2,row=3)
		ttk.Label(mainframe, text='The image dE_ds.png and the file dE_ds.txt have been generated').grid(column=2, row=4, sticky=W)
	except ValueError:
		print('ERROR IN FUNCTION : compute_plot_and_save_stopping_power')
		pass

def compute_Bragg_peak(material_var,z_var,m_var,Eps0,T_var):
# input  : material_var = name of the target material (ex. : "Boron (solid at room temperature)")
#          z_var        = projectile atomic number in ()
#          m_var        = projectile mass in (g)
#          Eps0         = initial projectile kinetic energy in (MeV)
#          T_var        = projectile stopping time in (s)
# output : [N,x,S]      = N  is the number of time iteration that is used to compute numerically
#                         dE/dt(t) = S[E(t)] v[E(t)], dx/dt(t) = v(E(t)) where (t_n = n dt with dt = T/N)
#                         x(t_n) in cm is the obtained array of successive positions and S[E(t_n)] the
#                         corresponding value of the stopping power in (MeV/cm)
	try :
		print('* Compute Bragg peak curve')
		N      = 1000 # compromise between accuracy and time computation
		dt     = (float(1.01*T_var) / float(N) )
		x      = np.zeros(N)
		v      = np.zeros(N)
		Eps    = np.zeros(N)
		S      = np.zeros(N)
		x[0]   = 0.
		Eps[0] = float(Eps0)
		v[0]   = physics_library.particle_velocity_norm(m_var,Eps[0])
		S[0]   = physics_library.Stopping_power(material_var,z_var,m_var,Eps[0])
		for k in range(0,N-1):
			Eps[k+1] = Eps[k] - ( S[k] * v[k] * dt )
			if (Eps[k+1] < 0.):
				break
			v[k+1]   = physics_library.particle_velocity_norm(m_var,Eps[k+1])
			S[k+1]   = max(0.,physics_library.Stopping_power(material_var,z_var,m_var,Eps[k+1]))
			x[k+1]   = x[k] + (v[k+1] * dt)
		return [N,x,S];
	except ValueError:
		print('ERROR IN FUNCTION : compute_Bragg_peak')
		pass

def plot_Bragg_peak(projectile_var,target_var,Eps0,R_var,x_var,S_var):
# input  : projectile_var = 'Proton' or 'Alpha particle'
#          target_var     = name of the target material (ex. : "Boron (solid at room temperature)")
#          Eps0           = initial projectile kinetic energy in (MeV)
#          R_var          = value of the projectile range in (cm)
#          x_var          = array containing all projectile successive position in (cm)
#          S_var          = array containing the corresponding value of the stopping power in (MeV/cm)
# output : nothing - it just plots the stopping power versus the 
#                    in flight projectile position 
#                and save the corresponding figure dE_dx.png
	try :
		# create image dE_dx.png
		font = {'style':  'normal',
		'color':  'black',
		'weight': 'normal',
		'size': 16,
				}
		xmax      = max(x_var) 
		n_xmax    = int(math.log(xmax)/math.log(10.)) - 1
		x_var     = x_var / 10.**n_xmax
		xmax_plot = float(int(12.5 * xmax / 10.**n_xmax)/10.)
		ymin      = S_var[0]
		for i in range(0,len(S_var)):
			if (S_var[i] <= ymin ) and (S_var[i] != 0.):
				ymin = S_var[i]
		n_ymin    = int(math.log(ymin)/math.log(10.)) - 1
		ymax      = max(S_var)
		n_ymax    = 1 + int(math.log(ymax)/math.log(10.))
		ymin_plot = float(10.**n_ymin)
		ymax_plot = float(10.**n_ymax)
		S2        = np.arange(ymin_plot,ymax_plot,(ymax_plot-ymin_plot)/500.)
		x2        = float(R_var / 10.**n_xmax) * np.ones(len(S2))
		fig=plt.figure()
		plt.semilogy(x_var,S_var ,'blue',linestyle='-' ,linewidth=2,label='E(x=0) = '+str(Eps0)+' MeV')
		plt.semilogy(x2   ,S2    ,'blue',linestyle='--',linewidth=2,label='R = '+"{0:.6f}".format(R_var)+' cm')
		plt.grid(True,which="both",axis="both",color='gray', linestyle='-', linewidth=1)
		plt.xticks(fontsize=16)
		plt.xlabel('x (1.e'+str(n_xmax)+' cm)', fontdict=font)
		plt.xlim([0,xmax_plot])
		plt.ylabel('S (MeV/cm)', fontdict=font)
		plt.yticks(fontsize=16)
		plt.ylim([ymin_plot,ymax_plot])
		ttl = str(projectile_var)+' Bragg peak in \n'+str(target_var)
		plt.title(ttl, fontdict=font)
		leg = plt.legend(loc='upper left', fontsize=16, fancybox=True)
		leg.get_frame().set_alpha(0.5)
		fig.savefig("dE_dx.png",bbox_inches='tight')
		print('* Bragg peak plot is saved')
	except ValueError:
		print('ERROR IN FUNCTION : plot_Bragg_peak')
		print('ValueError = '+str(ValueError))
		pass

def save_Bragg_peak(N_gridpoints_var,x_var,S_var) :
# input  : N_gridpoints_var = Number of gridpoints
#          x_var            = array containing all projectile successive positions in (cm)
#          S_var            = array containing the corresponding value of the stopping power in (MeV/cm)
# output : nothing - it just saves the text file dE_dx.txt containing the stopping power versus the 
#                    projectile position
	try :
		# create dE_dx.txt file
		txtfile = open("dE_dx.txt","w+")
		txtfile.write('      x (cm)     '+' '+'   S (MeV/cm)    '+" \n")
		for i in range(0,N_gridpoints_var):
			x_str = "{0:.15f}".format(x_var[i]) 
			S_str   = "{0:.15f}".format(S_var[i])
			txtfile.write(x_str+' '+S_str+" \n")
		txtfile.close()
		print('* Bragg peak text file is saved')
		print('-------------------------------')
	except ValueError:
		print('ERROR IN FUNCTION : save_Bragg_peak')
		pass

def compute_range(projectile_var,target_var,Eps_ini_var):
# input  : projectile_var = 'Proton' or 'Alpha particle'
#          target_var     = name of the target material (ex. : "Boron (solid at room temperature)")
#          Eps_ini_var    = value of projectile initial kinetic energy in (MeV)
# output : [R,T]          = Projectile range R in (cm) and corresponding stopping time T in (S)
	try :
		print('* Compute projectile range')
		# projectile
		if ( str(projectile_var) == 'Proton' ) :
			z = 1
			m = mp
		elif ( str(projectile_var) == 'Alpha particle' ) :
			z = 2
			m = ma
		# target
		material = get_material(target_var)
		# Compute the projectile range R and Stopping time T
		dEps         = float(Eps_ini_var)*1.e-3
		N_gridpoints = 1 + int(float(Eps_ini_var)/dEps)
		dEps_var     = float(Eps_ini_var) / N_gridpoints
		Eps          = np.zeros(N_gridpoints)
		x            = np.zeros(N_gridpoints)
		v            = np.zeros(N_gridpoints)
		S            = np.zeros(N_gridpoints)
		R            = 0.
		T            = 0.
		for i in range(1,N_gridpoints):
			Eps[i] = float(i)*float(dEps_var)
			S[i]   = physics_library.Stopping_power(material,z,m,Eps[i])
			v[i]   = physics_library.particle_velocity_norm(m,Eps[i])
			T      = T + ( float(dEps_var) / ( v[i] * S[i]) )
			R      = R + ( float(dEps_var) / S[i] )
		return [R,T];
	except ValueError:
		print('ERROR IN FUNCTION : compute_range')
		pass

def compute_plot_and_save_Bragg_peak(projectile_var,target_var,Eps_ini_var):
# input  : projectile_var = 'Proton' or 'Alpha particle'
#          target_var     = name of the target material (ex. : "Boron (solid at room temperature)")
#          Eps_ini_var    = projectile initial kinetic energy in (MeV)
# output : nothing - it just plots the stopping power versus the 
#                    in flight projectile position, save the corresponding figure dE_dx.png, 
#                    and generate the corresponding text file dE_dx.txt containing all values 
	try:
		# projectile
		if ( str(projectile_var) == 'Proton' ) :
			z = 1
			m = mp
		elif ( str(projectile_var) == 'Alpha particle' ) :
			z = 2
			m = ma
		# target
		material = get_material(target_var)
		# Compute the projectile range R and Stopping time T
		[R,T] = compute_range(projectile_var,target_var,Eps_ini_var)
		ttk.Label(mainframe, text='R = '+str(R)+' cm and Stopping Time : T = '+str(T)+ ' s').grid(column=4, row=2, sticky=W)
		# compute Bragg peak
		nrj_0 = float(Eps_ini_var)
		res   = compute_Bragg_peak(material,z,m,nrj_0,T)
		N     = res[0]
		x     = res[1]
		S     = res[2]
		# plot
		plot_Bragg_peak(projectile_var,target_var,nrj_0,R,x,S)
		# save
		save_Bragg_peak(N,x,S)
		# display
		img = Image.open('dE_dx.png')
		img = img.resize((460, 355), Image.ANTIALIAS)
		img = ImageTk.PhotoImage(img)
		ttk.Label(mainframe, image = img).grid(column=4,row=3)
		ttk.Label(mainframe, text='The image dE_dx.png and the file dE_dx.txt have been generated').grid(column=4, row=4, sticky=E)
	except ValueError:
		print('ERROR IN FUNCTION : compute_plot_and_save_Bragg_peak')
		pass

###################

root = Tk()

# Frame's title
root.title("Stopping Powers of Protons and Alpha particles in Ambient Matter (SPAM) by Michael J TOUATI - 09/27/2019 - https://github.com/michaeltouati")
# Frame's properties
mainframe = ttk.Frame(root, padding="200 200 200 200")
mainframe.grid(column=2, row=11, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

# Projectile type entry
ttk.Label(mainframe, text="Projectile :").grid(column=1, row=1, sticky=E)
projectile = StringVar()
projectile_entry = ttk.Combobox(mainframe, width=12, values=['Proton', 'Alpha particle'])
projectile_entry.grid(column=2, row=1, sticky=(W, E))

# Target type entry
ttk.Label(mainframe, text="Target :").grid(column=1, row=2, sticky=E)
target       = StringVar()
table        = generate_material_list()
target_entry = ttk.Combobox(mainframe, width=50, values=table)
target_entry.grid(column=2, row=2, sticky=(W, E))

# Left button
Button_left_entry = ttk.Button(mainframe, text="Plot / Save", command=calculate_left).grid(column=1, row=3, sticky=E)

# Entry of intial projectile kinetic energy in (MeV)
ttk.Label(mainframe, text="Kinetic Energy in (MeV) :").grid(column=3, row=1, sticky=E)
Eps_ini = StringVar()
Eps_ini_entry = ttk.Entry(mainframe, width=50, textvariable=Eps_ini)
Eps_ini_entry.grid(column=4, row=1, sticky=(W, E))

# Deduced range
ttk.Label(mainframe, text="Deduced Range :").grid(column=3, row=2, sticky=E)

# Righ button
Button_right_entry = ttk.Button(mainframe, text="Plot / Save", command=calculate_right).grid(column=3, row=3, sticky=E)

# Quit button
Button_quit = ttk.Button(mainframe, text="Quit", command=close).grid(column=3, row=4, sticky=(W,E))

# Frame reactualization
for child in mainframe.winfo_children(): 
	child.grid_configure(padx=2, pady=11)
projectile_entry.focus()
target_entry.focus()
Eps_ini_entry.focus()
root.mainloop()
