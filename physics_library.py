#######################################################################
##                                                                   ##
##  Stopping Power of Protons and Alpha particles in Ambient Matter  ##
##                              (SPAM)                               ##
##                                                                   ##
## Copyright © 2020 Michaël J TOUATI                                 ##
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
import math
import maths_library
import physical_constants

pi    = math.pi
mu    = physical_constants.atomic_mass_unit
mp    = physical_constants.proton_mass
me    = physical_constants.electron_mass     
qe    = physical_constants.elementary_charge
c     = physical_constants.light_speed_in_vacuum
mpc2  = physical_constants.proton_mass_energy
kB    = physical_constants.Boltzmann_constant
h     = physical_constants.Planck_constant
hbar  = h / (2.*pi)
alpha = physical_constants.fine_structure_constant
eV    = physical_constants.eV_in_erg
keV   = physical_constants.keV_in_erg
MeV   = physical_constants.MeV_in_erg
Ry    = physical_constants.Rydberg_constant
s     = physical_constants.symbol

def get_atomic_number(material):
# input  : material = material name (string of characters)
# output : Z        = material atomic number in ()
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				Z = n + 1
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					Z = n + 1
					found = True
	if (found == False):
		if (material == 'GafChromic HD-V2 (solid at room temperature)') :
			Z = 7.58 # (effective charge)
		elif (material == 'Polyester (solid at room temperature)'):
			Z = 6.64
		elif (material == 'Plastic scintillator (vinyltoluene-based solid at room temperature)'):
			Z = 3.383
		elif (material == 'Water (liquid at room temperature)') :
			Z = 7.216742
		elif (material == 'Air (dry gas at room temperature)') :
			Z = 7.372747
		elif (material == 'Glass SiO2 (solid at room temperature)') :
			Z = 10.804602
		else :
			Z = 0
	return Z;

def get_density(material):
# input  : material = material name (string of characters)
# output : rho      = material density in (g/cm^3)
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				rho = eval("physical_constants.density_"+s[n])
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					rho = eval("physical_constants.density_"+s[n])[m]
					found = True
	if (found == False):
		if (material == 'GafChromic HD-V2 (solid at room temperature)'):
			rho = 1.2
		elif (material == 'Polyester (solid at room temperature)'):
			rho = 1.35
		elif (material == 'Plastic scintillator (vinyltoluene-based solid at room temperature)'):
			rho = 1.032
		elif (material == 'Water (liquid at room temperature)'):
			rho = 1.
		elif (material == 'Air (dry gas at room temperature)') :
			rho = 1.20484e-3
		elif (material == 'Glass SiO2 (solid at room temperature)') :
			rho = 2.32
		else :
			rho = 0
	return rho;

def get_mean_excitation_energy(material):
# input  : material = material name (string of characters)
# output : I        = mean excitation energy in (eV)
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				I = eval("physical_constants.mean_excitation_energy_"+s[n])
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					I = eval("physical_constants.mean_excitation_energy_"+s[n])[m]
					found = True
	if (found == False):
		print("This is a compound")
		I = 0
	return I;

def get_standard_atomic_weight(material):
# input  : material = material name (string of characters)
# output : A        = material standard atomic weight in (g/mol)
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				A = eval("physical_constants.standard_atomic_weight_"+s[n])
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					A = eval("physical_constants.standard_atomic_weight_"+s[n])
					found = True
	if (found == False):
		print("This is a compound")
		A = 0
	return A;

def get_low_energy_proton_empirical_formula_A_parameters(material):
# input  : material = material name (string of characters)
# output : A        = material standard atomic weight in (g/mol)
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				A1  = eval("physical_constants.coeff_A1_low_energy_proton_empirical_formula_"+s[n])
				A2  = eval("physical_constants.coeff_A2_low_energy_proton_empirical_formula_"+s[n])
				A3  = eval("physical_constants.coeff_A3_low_energy_proton_empirical_formula_"+s[n])
				A4  = eval("physical_constants.coeff_A4_low_energy_proton_empirical_formula_"+s[n])
				A5  = eval("physical_constants.coeff_A5_low_energy_proton_empirical_formula_"+s[n])
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					A1  = eval("physical_constants.coeff_A1_low_energy_proton_empirical_formula_"+s[n])[m]
					A2  = eval("physical_constants.coeff_A2_low_energy_proton_empirical_formula_"+s[n])[m]
					A3  = eval("physical_constants.coeff_A3_low_energy_proton_empirical_formula_"+s[n])[m]
					A4  = eval("physical_constants.coeff_A4_low_energy_proton_empirical_formula_"+s[n])[m]
					A5  = eval("physical_constants.coeff_A5_low_energy_proton_empirical_formula_"+s[n])[m]
					found = True
	if (found == False):
		print("This is a compound")
		A1 = 0
		A2 = 0
		A3 = 0
		A4 = 0
		A5 = 0
	return [A1,A2,A3,A4,A5];

def get_low_energy_proton_empirical_formula_T_parameters(material):
# input  : material = material name (string of characters)
# output : A        = material standard atomic weight in (g/mol)
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				T1  = eval("physical_constants.param_T1_low_energy_proton_empirical_formula_"+s[n])
				T2  = eval("physical_constants.param_T2_low_energy_proton_empirical_formula_"+s[n])
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					T1  = eval("physical_constants.param_T1_low_energy_proton_empirical_formula_"+s[n])[m]
					T2  = eval("physical_constants.param_T2_low_energy_proton_empirical_formula_"+s[n])[m]
					found = True
	if (found == False):
		print("This is a compound")
		T1 = 0
		T2 = 0
	return [T1,T2];

def get_low_energy_alpha_empirical_formula_a_parameters(material):
# input  : material = material name (string of characters)
# output : A        = material standard atomic weight in (g/mol)
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				a1  = eval("physical_constants.coeff_a1_low_energy_alpha_empirical_formula_"+s[n])
				a2  = eval("physical_constants.coeff_a2_low_energy_alpha_empirical_formula_"+s[n])
				a3  = eval("physical_constants.coeff_a3_low_energy_alpha_empirical_formula_"+s[n])
				a4  = eval("physical_constants.coeff_a4_low_energy_alpha_empirical_formula_"+s[n])
				a5  = eval("physical_constants.coeff_a5_low_energy_alpha_empirical_formula_"+s[n])
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					a1  = eval("physical_constants.coeff_a1_low_energy_alpha_empirical_formula_"+s[n])[m]
					a2  = eval("physical_constants.coeff_a2_low_energy_alpha_empirical_formula_"+s[n])[m]
					a3  = eval("physical_constants.coeff_a3_low_energy_alpha_empirical_formula_"+s[n])[m]
					a4  = eval("physical_constants.coeff_a4_low_energy_alpha_empirical_formula_"+s[n])[m]
					a5  = eval("physical_constants.coeff_a5_low_energy_alpha_empirical_formula_"+s[n])[m]
					found = True
	if (found == False):
		print("This is a compound")
		a1 = 0
		a2 = 0
		a3 = 0
		a4 = 0
		a5 = 0
	return [a1,a2,a3,a4,a5];

def get_low_energy_alpha_empirical_formula_T_parameters(material):
# input  : material = material name (string of characters)
# output : A        = material standard atomic weight in (g/mol)
	N     = len(s)
	found = False
	for n in range(0,N) :
		mat = eval("physical_constants.name_"+s[n])
		if ( len(mat) > 10 ) : # no elements with more than 10 different phases
			if ( mat == material ):
				T1  = eval("physical_constants.param_T1_low_energy_alpha_empirical_formula_"+s[n])
				T2  = eval("physical_constants.param_T2_low_energy_alpha_empirical_formula_"+s[n])
				found = True
		else :
			for m in range(0,len(mat)):
				mat2 = mat[m]
				if ( mat2 == material ) :
					T1  = eval("physical_constants.param_T1_low_energy_alpha_empirical_formula_"+s[n])[m]
					T2  = eval("physical_constants.param_T2_low_energy_alpha_empirical_formula_"+s[n])[m]
					found = True
	if (found == False):
		print("This is a compound")
		T1 = 0
		T2 = 0
	return [T1,T2];

def Stopping_power(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
	if (material == 'GafChromic HD-V2 (solid at room temperature)'):
		S = Stopping_power_in_compounds(material,z,m,nrj)
	elif (material == 'Polyester (solid at room temperature)'):
		S = Stopping_power_in_compounds(material,z,m,nrj)
	elif (material == 'Plastic scintillator (vinyltoluene-based solid at room temperature)'):
		S = Stopping_power_in_compounds(material,z,m,nrj)
	elif (material == 'Water (liquid at room temperature)') :
		S = Stopping_power_in_compounds(material,z,m,nrj)
	elif (material == 'Air (dry gas at room temperature)') :
		S = Stopping_power_in_compounds(material,z,m,nrj)
	elif (material == 'Glass SiO2 (solid at room temperature)') :
		S = Stopping_power_in_compounds(material,z,m,nrj)
	else :
		if (int(z) == 1) :
			[T1,T2] = get_low_energy_proton_empirical_formula_T_parameters(material)
			Z_mat = get_atomic_number(material)
			if (Z_mat != 29) and (Z_mat != 32) and (Z_mat != 36) and (Z_mat != 47) and (Z_mat != 51) and (Z_mat != 79):
				T1 = T2 
			if (nrj <= T1):
				S = Stopping_power_low_energy_proton(material, nrj)
			elif (nrj > T2):
				S = Stopping_power_Bethe(material,z,m,nrj)
			else :
				S1    = Stopping_power_low_energy_proton(material, T1)
				S2    = Stopping_power_Bethe(material,z,m,T2)
				S     = maths_library.linear_interpolation_1D(T1,S1,T2,S2,nrj)
			S = S + Stopping_power_nuclear_Moliere(material,z,m,nrj)
		elif (int(z) == 2) :
			[T1,T2] = get_low_energy_alpha_empirical_formula_T_parameters(material)
			if (nrj <= T1):
				S = Stopping_power_low_energy_alpha(material, nrj)
			elif (nrj > T2):
				S = Stopping_power_Bethe(material,z,m,nrj)
			else :
				beta  = particle_beta_norm(m,nrj)
				beta2 = beta * beta
				S1    = beta2 * Stopping_power_low_energy_alpha(material, T1)
				S2    = beta2 * Stopping_power_Bethe(material,z,m,T2)
				S     = maths_library.linear_interpolation_1D(T1,S1,T2,S2,nrj)/beta2
			S = S + Stopping_power_nuclear_Ziegler(material,z,m,nrj)
		else :
			# if (nrj <= 1.e1) :
			# 	S = Stopping_power_low_energy(Z,z,m,nrj)
			# else :
			# 	S = Stopping_power_Bethe(Z,z,m,nrj)
			S = Stopping_power_Bethe(material,z,m,nrj)
			S = S + Stopping_power_nuclear_Ziegler(material,z,m,nrj) 
	return S;

######################################################################################
#                                 Special relativity                                 #
######################################################################################

def particle_Lorentz_factor(m,nrj):
# input  : m     = particle mass m in g
#          nrj   = particle kinetic energy in (MeV)
# output : g     = particle Lorentz factor in ()
	mc2 = m * c * c
	E   = nrj * MeV
	g   = 1. + ( E / mc2 )
	return g;

def particle_beta_norm(m,nrj):
# input  : m    = particle mass m in g
#          nrj  = particle kinetic energy in (MeV)
# output : beta = v / c in () where v is the particle velocity
	g = particle_Lorentz_factor(m,nrj)
	b = math.sqrt( 1. - ( 1. / (g**2.) ) )
	return b;

def particle_velocity_norm(m,nrj):
# input  : m   = particle mass m in g
#          nrj = particle kinetic energy in (MeV)
# output : v   = particle velocity in (cm/s)
	b = particle_beta_norm(m,nrj)
	v = c * b
	return v;

def particle_momentum_norm(m,nrj):
# input  : m   = particle mass m in g
#          nrj = particle kinetic energy in (MeV)
# output : p   = particle momentum in (g.cm/s)     
	g = particle_Lorentz_factor(m,nrj)
	v = particle_velocity_norm(m,nrj)
	p = g * m * v
	return p;

######################################################################################
#                                 Bethe Stopping number                              #
# _ Bethe, H. A. - Quantenmechanik der Ein-und zei-elektronen Preblem -              #
#   Handbuch der Physik Vol. 24 , No. 1, p. 273                                      #
# _ Bethe H. A. and Ashkin J. - Passage of Radiations through matter -               # 
#   Exp. Nucl. Phys., Vol. 1, Segré E. Ed.                                           #
######################################################################################

def largest_possible_energy_loss_in_single_collision(m,nrj):
# input  : m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Wm = largest possible energy loss in (erg) :
	mec2  = me * c * c
	beta  = particle_beta_norm(m,nrj)
	if (beta > 0.01):
		beta2 = beta * beta
		y     = me / m
		x     = 1. - beta2
		D     = 1. + ( 2. * y / (x**0.5) ) + (y**2.)
		Wm    = 2. * mec2 * beta2 / ( x * D )
	else :
		v  = c * beta
		v2 = v * v
		Wm = 2. * me * v2 
	return Wm;

def Stopping_number_Bethe(material,m,nrj):
# input  : material = material name (string of characters)
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : Bethe stopping number L0(beta) in () :
	mec2  = me * c * c
	beta  = particle_beta_norm(m,nrj)
	beta2 = beta * beta
	x     = 1. - beta2
	Wm    = largest_possible_energy_loss_in_single_collision(m,nrj)
	I     = get_mean_excitation_energy(material) * eV
	L0    = ( 0.5 * math.log( 2. * mec2 * Wm *beta2 / x ) ) - math.log(I) - beta2
	return L0;

def Stopping_power_Bethe(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
	# Bethe stopping number
	L0 = Stopping_number_Bethe(material,m,nrj)
	Z  = get_atomic_number(material)
	# Shell correction
	C  = Stopping_number_shell_correction(Z,m,nrj)
	# We neglect density effects assuming non relativistic projectile energy
	D = 0.
	L0 = L0 - C - D
	# Barkas correction
	L1 = Stopping_number_Barkas_correction(Z,z,m,nrj)
	# Bloch correction
	L2 = Stopping_number_Bloch_correction(z,m,nrj)
	#
	L  = L0 + L1 + L2
	# electron density
	rho = get_density(material)
	A   = get_standard_atomic_weight(material)
	mi  = A * mu
	ni  = rho / mi
	ne  = float(Z) * ni
	# factor
	v     = particle_velocity_norm(m, nrj)
	v2    = v * v
	qe2   = qe * qe
	qe4   = qe2 * qe2
	z2    = float(z) * float(z)
	alpha = 4. * pi * ne * qe4 * z2 / (me * v2)
	#
	S = alpha * L / MeV
	return S;

######################################################################################
#                                    Density effect                                  #
#                        W. Swann J. Franklin lnst. 226, 598 (1938)                  #
#                            E. Fermi Phys. Rev. 57,485 (1940)                       #
######################################################################################

def delta_correction(material,m,nrj):
# input  : material = material name (string of characters)
#          m     = projectile mass in g
#          nrj   = projectile energy in (MeV)
# output : delta = Stopping number correction such that L0(beta) -> L0(beta) - (delta/2) 
#          valid only when beta -> 1
	# plasma frequency
	rho = get_density(material)
	A   = get_standard_atomic_weight(material)
	mi  = A * mu
	ni  = rho / mi
	ne  = float(Z) * ni
	omega_pe = math.sqrt(4. * pi * ne * (qe**2.) / me)
	# mean excitation potential
	I = get_mean_excitation_energy(material) * eV
	#
	alpha  = (hbar * omega_pe / I)
	alpha2 = alpha * alpha
	beta   = particle_beta_norm(m,nrj)
	beta2  = beta * beta
	delta  = math.log( alpha2 / (1.-beta2) ) - 1.
	return delta;

######################################################################################
#                                 Barkas correction                                  #
#     Barkas et al. - Resolution of the Σ- mass anomaly - PRL Vol. II, No. 1, p. 6   #
######################################################################################

def Barkas_correction(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction : 
#          L1(beta) = g F_ARB( b/sqrt(x) ) / ( sqrt(Z) * x^(3/2) ) in ()
#          where g = constant = 1.29 and x = (beta/alpha)^2 / Z with alpha the 
#          fine structure constant according to Eq (2.8) of Berger et al. 
#          - Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
	g     = 1.29
	beta  = particle_beta_norm(m,nrj)
	s     = physical_constants.symbol
	b     = eval("physical_constants.scaled_minimum_impact_parameter_"+s[int(Z)-1])
	x     = ( ( beta / alpha )**2. ) / float(Z)
	w = b / (x**0.5)
	F_ARB = maths_library.Ashley_function(w)
	L1    = g * F_ARB / ( ( Z * (x**3.) )**0.5 )
	return L1;

def Bichsel_correction(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction for int(Z) >63: 
#          L1(beta) = g1 beta^(-2*g2) in ()
#          according to H. Bichsel - Phys. Rev. A, Vol. 41, No. 7, p. 3642 (1990)
	g1    = 0.002833
	g2    = 0.6
	beta  = particle_beta_norm(m,nrj)
	L1    = g1 / ( beta**(2.*g2) )
	return L1;

def Bichsel_correction_silver(m,nrj):
# input  : m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction for Silver : 
#          L1(beta) = g1 beta^(-2*g2) in ()
#          according to H. Bichsel - Phys. Rev. A, Vol. 41, No. 7, p. 3642 (1990)
	g1    = 0.006812
	g2    = 0.45
	beta  = particle_beta_norm(m,nrj)
	L1    = g1 / ( beta**(2.*g2) ) 
	return L1;


def Stopping_number_Barkas_correction(Z,z,m,nrj):
# input  : Z   = Material atomic number in ()
#          z   = projectile atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction : L1(beta) in ()
#          according to Berger et al. - 
#          - Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993) 
	if (int(Z) < 64) and int(Z) != 47:
		L1 = Barkas_correction(Z,m,nrj)
	elif int(Z) == 47:
		L1 = Bichsel_correction_silver(m,nrj)
	else :
		L1 = Bichsel_correction(Z,m,nrj)
	L1 = float(z) * L1 
	return L1;

######################################################################################
#                                 Bloch correction                                   #
#       F. Bloch - Zur Bremsung rash bewegter teilchen beim Durchgang durch die      #
#                         Matter - Ann. Phys. 16, 285 (1933)                         #
######################################################################################

def Stopping_number_Bloch_correction(z,m,nrj):
# input  : z   = projectile atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction 
#                           ----infinity
#                           \                   1
#          L2(beta) = -y^2 * |           --------------- in ()
#                           /             n ( n^2 + y^2)
#                           ----n=1
#          according to Berger et al. - 
#          Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
	beta  = particle_beta_norm(m,nrj)
	y     = float(z) * alpha / beta
	L2    = maths_library.Bloch_function(y)
	return L2;

######################################################################################
#                                 Shell corrections                                  #
#      M. S. Livingston and H. A. Bethe, Revs. Modern Phys. 9, 263 (1937)            #
######################################################################################

def theta_K(Z):
# input  : Z     = Atomic number in ()
# output : theta = 1 - (I_R - I_0 )/I_NR where I_NR = theoretical non-relativistic ionization potential of K-shell electron(s)
#                                              I_R  = theoretical relativistic ionization potential of K-shell electron(s)
#                                              I_0  = experimental energy difference from the K-shell to the lowest-lying unoccupied state
#   M. C. Walske, Phys. Rev., Vol. 101, No. 3, p. 940 (1956) and
#   V. H. Hönl, Z. Physik 84, 1 (1933) :
	theta_K_array = [0.        , 0.        , 0.51261242, 0.57623538, 0.61299021, 0.63698434, 0.65390287, 0.66650395, 0.68190   , 0.69100   , 0.70010   , 0.70920   , 0.71825   , 0.72619   , 0.73359   , 0.74044   ,    0.74626, 0.75143   , 0.75740   , 0.76229   , 0.76637   , 0.77030   , 0.77424   , 0.77779   , 0.78155   , 0.78474   , 0.78791   , 0.79094   , 0.79375   , 0.79656   , 0.79874   , 0.80043   , 0.80358   , 0.80667   , 0.80781   , 0.80968   , 0.81255   , 0.81431   , 0.81566   , 0.81735   , 0.81903   , 0.82071   , 0.82239   , 0.82408   , 0.82559   , 0.82703   , 0.82884   , 0.82992   , 0.83213   , 0.83278   , 0.83556   , 0.83591   , 0.83800   , 0.83905   , 0.83978   , 0.84189   , 0.84248   , 0.84360   , 0.84472   , 0.84584   , 0.84639   , 0.84860   , 0.84922   , 0.84925   , 0.85132   , 0.85255   , 0.85267   , 0.85266   , 0.85361   , 0.85519   , 0.85597   , 0.85608   , 0.85694   , 0.85890  , 0.85950   , 0.85949   , 0.85948   , 0.86099   , 0.86250   , 0.86290   , 0.86290   , 0.86289]
	z = int(Z) - 1
	theta = theta_K_array[z]
	return theta;

def S_K(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_K(Z,eta_K) = S(Z) log(eta_K) + T(Z) - C_K(Z,eta_K)
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	S_K_array = [2.04646, 2.01633, 1.98708, 1.95867, 1.93424, 1.91236, 1.89134, 1.87105, 1.85107, 1.83228, 1.81350, 1.79608, 1.77875, 1.76431, 1.75107, 1.73883, 1.72877, 1.72001, 1.71046, 1.70273, 1.69639, 1.69028, 1.68415, 1.67863, 1.67294, 1.66830, 1.66369, 1.65928, 1.65519, 1.65110, 1.64793, 1.64551, 1.64123, 1.63702, 1.63547, 1.63293, 1.62903, 1.62663, 1.62480, 1.62250, 1.62021, 1.61799, 1.61584, 1.61367, 1.61174, 1.60990, 1.60758, 1.60620, 1.60337, 1.60254, 1.59898, 1.59853, 1.59586, 1.59451, 1.59358, 1.59099, 1.59027, 1.58890, 1.58754, 1.58617, 1.58550, 1.58280, 1.58205, 1.58201, 1.57954, 1.57809, 1.57794, 1.57796, 1.57684, 1.57497, 1.57405, 1.57392, 1.57291, 1.57059, 1.56988, 1.56990, 1.56991, 1.56818, 1.56647, 1.56602, 1.56602, 1.56603]
	z         = int(Z) - 1
	S         = S_K_array[z]
	return S; 

def T_K(Z):
# input  : Z = Atomic number in ()
# output : T = T coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_K(Z,eta_K) = S(Z) log(eta_K) + T(Z) - C_K(Z,eta_K)
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	T_K_array = [2.55340, 2.54430, 2.53520, 2.52610, 2.51687, 2.50795, 2.49880, 2.48948, 2.48007, 2.47038, 2.46069, 2.45068, 2.44072, 2.43183, 2.42351, 2.41577, 2.40885, 2.40281, 2.39618, 2.39064, 2.38591, 2.38135, 2.37678, 2.37266, 2.36827, 2.36453, 2.36080, 2.35724, 2.35394, 2.35064, 2.34808, 2.34609, 2.34235, 2.33869, 2.33734, 2.33512, 2.33172, 2.32964, 2.32804, 2.32604, 2.32404, 2.32205, 2.32006, 2.31806, 2.31627, 2.31456, 2.31242, 2.31114, 2.30852, 2.30775, 2.30446, 2.30404, 2.30157, 2.30032, 2.29946, 2.29693, 2.29622, 2.29488, 2.29353, 2.29219, 2.29153, 2.28887, 2.28813, 2.28810, 2.28562, 2.28416, 2.28402, 2.28403, 2.28290, 2.28102, 2.28009, 2.27996, 2.27894, 2.27660, 2.27589, 2.27590, 2.27591, 2.27412, 2.27232, 2.27184, 2.27184, 2.27186]
	z         = int(Z) - 1
	T         = T_K_array[z]
	return T; 

def U_K(Z):
# input  : Z = Atomic number in ()
# output : U = U coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              C_K(Z,eta_K) = (U_K / eta_K) + (V_K / eta_K^2)
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	U_K_array = [1.96137, 1.97275, 1.98412, 1.99549, 2.00731, 2.01909, 2.02986, 2.03987, 2.04954, 2.05792, 2.06627, 2.07333, 2.08034, 2.08566, 2.09039, 2.09473, 2.09788, 2.10061, 2.10360, 2.10585, 2.10754, 2.10917, 2.11080, 2.11228, 2.11370, 2.11474, 2.11577, 2.11675, 2.11766, 2.11858, 2.11929, 2.11980, 2.12057, 2.12133, 2.12161, 2.12207, 2.12277, 2.12320, 2.12353, 2.12395, 2.12436, 2.12472, 2.12500, 2.12529, 2.12555, 2.12579, 2.12610, 2.12628, 2.12666, 2.12677, 2.12724, 2.12730, 2.12766, 2.12783, 2.12796, 2.12822, 2.12829, 2.12843, 2.12856, 2.12870, 2.12876, 2.12903, 2.12910, 2.12910, 2.12931, 2.12942, 2.12944, 2.12943, 2.12952, 2.12966, 2.12973, 2.12974, 2.12982, 2.13000, 2.13005, 2.13005, 2.13005, 2.13014, 2.13021, 2.13023, 2.13023, 2.13023]
	z         = int(Z) - 1
	U         = U_K_array[z]
	return U; 

def V_K(Z):
# input  : Z = Atomic number in ()
# output : V = V coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              C_K(Z,eta_K) = (U_K / eta_K) + (V_K / eta_K^2)
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	V_K_array = [8.35644, 8.35189, 8.34734, 8.34280, 8.33896, 8.33578, 8.33301, 8.33060, 8.32832, 8.32649, 8.32468, 8.32341, 8.32214, 8.32134, 8.32067, 8.32006, 8.31966, 8.31935, 8.31917, 8.31906, 8.31900, 8.31894, 8.31888, 8.31883, 8.31882, 8.31887, 8.31891, 8.31896, 8.31900, 8.31904, 8.31908, 8.31911, 8.31924, 8.31936, 8.31941, 8.31948, 8.31960, 8.31967, 8.31972, 8.31979, 8.31986, 8.31994, 8.32004, 8.32014, 8.32023, 8.32032, 8.32043, 8.32049, 8.32062, 8.32066, 8.32083, 8.32085, 8.32097, 8.32104, 8.32108, 8.32123, 8.32127, 8.32135, 8.32143, 8.32150, 8.32154, 8.32170, 8.32174, 8.32174, 8.32190, 8.32200, 8.32201, 8.32201, 8.32208, 8.32221, 8.32227, 8.32228, 8.32235, 8.32251, 8.32255, 8.32255, 8.32255, 8.32268, 8.32282, 8.32286, 8.32286, 8.32286]
	z         = int(Z) - 1
	V         = V_K_array[z]
	return V; 

def B_K(theta,eta):
# input  : theta = 1 - (I_R - I_0 )/I_NR where I_NR = theoretical non-relativistic ionization potential of L-shell electron(s)
#                                              I_R  = theoretical relativistic ionization potential of L-shell electron(s)
#                                              I_0  = experimental energy difference from the L-shell to the lowest-lying unoccupied state
#          eta   = me v^2 / (2 Zeff^2 Ry) where v is the projectile velocity and Zeff the effective charge seen by a K-shell electron
# output : B     = K-shell contribution to the stopping number
#
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 : B_K_table[i][j] = B_K(eta_K[i],theta_K[j]) with :
	theta_K_array = [0.95, 0.94, 0.92, 0.9, 0.88, 0.86, 0.85, 0.84, 0.82, 0.8, 0.78, 0.76, 0.75, 0.74, 0.72, 0.7, 0.68, 0.66, 0.65, 0.64]
	eta_K_array   = [0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.5, 1.7, 2.0, 3.0, 5.0, 7.0, 10.0]
	B_K_table     = [[1.34782e-08, 1.46132e-08, 1.72179e-08, 2.03521e-08, 2.41370e-08, 2.87247e-08, 3.13778e-08, 3.43072e-08, 4.11274e-08, 4.94946e-08, 5.98040e-08, 7.25636e-08, 8.00602e-08, 8.84294e-08, 1.08253e-07, 1.33148e-07, 1.64573e-07, 2.04459e-07, 2.28346e-07, 2.55370e-07], [6.87555e-08, 7.44373e-08, 8.74397e-08, 1.03022e-07, 1.21760e-07, 1.44370e-07, 1.57398e-07, 1.71747e-07, 2.05023e-07, 2.45620e-07, 2.95345e-07, 3.56497e-07, 3.92247e-07, 4.32017e-07, 5.25688e-07, 6.42391e-07, 7.88464e-07, 9.72171e-07, 1.08140e-06, 1.20435e-06], [3.78413e-07, 4.08831e-07, 4.78154e-07, 5.60760e-07, 6.59478e-07, 7.77847e-07, 8.45709e-07, 9.20187e-07, 1.09192e-06, 1.29981e-06, 1.55232e-06, 1.86011e-06, 2.03881e-06, 2.23662e-06, 2.69889e-06, 3.26860e-06, 3.26860e-06, 4.84882e-06, 5.36428e-06, 5.94048e-06], [2.53200e-06, 2.72664e-06, 3.16738e-06, 3.68791e-06, 4.30423e-06, 5.03578e-06, 5.45200e-06, 5.90633e-06, 6.94501e-06, 8.18757e-06, 9.67802e-06, 1.14707e-05, 1.25008e-05, 1.36329e-05, 1.62480e-05, 1.94200e-05, 2.32783e-05, 2.79850e-05, 3.07181e-05, 3.37432e-05], [9.43891e-06, 1.01339e-05, 1.16984e-05, 1.35313e-05, 1.56829e-05, 1.82138e-05, 1.96439e-05, 2.11973e-05, 2.47216e-05, 2.88935e-05, 3.38425e-05, 3.97259e-05, 4.30763e-05, 4.67351e-05, 5.51033e-05, 6.51154e-05, 7.71154e-05, 9.15431e-05, 9.98212e-05, 1.08909e-04], [5.67088e-05, 6.05576e-05, 6.91311e-05, 7.90324e-05, 9.04832e-05, 1.03744e-04, 1.11147e-04, 1.19122e-04, 1.36980e-04, 1.57744e-04, 1.81920e-04, 2.10106e-04, 2.25918e-04, 2.43007e-04, 2.81460e-04, 3.26458e-04, 3.79175e-04, 4.41006e-04, 4.75845e-04, 5.13606e-04], [1.91576e-04, 2.03626e-04, 2.30230e-04, 2.60584e-04, 2.95248e-04, 3.34870e-04, 3.56771e-04, 3.80200e-04, 4.32104e-04, 4.91584e-04, 5.59802e-04, 6.38103e-04, 6.81511e-04, 7.28042e-04, 8.31425e-04, 9.50341e-04, 1.08721e-03, 1.24485e-03, 1.33245e-03, 1.42650e-03], [4.74226e-04, 5.02006e-04, 5.62872e-04, 6.31597e-04, 7.09244e-04, 7.97023e-04, 8.45134e-04, 8.96410e-04, 1.00867e-03, 1.13590e-03, 1.28002e-03, 1.44336e-03, 1.53305e-03, 1.62855e-03, 1.83861e-03, 2.07396e-03, 2.34750e-03, 2.65469e-03, 2.82358e-03, 3.00358e-03], [9.67285e-04, 1.02029e-03, 1.13566e-03, 1.26476e-03, 1.46928e-03, 1.57113e-03, 1.65921e-03, 1.75244e-03, 1.95562e-03, 2.18336e-03, 2.42872e-03, 2.72510e-03, 2.88111e-03, 3.04636e-03, 3.40681e-03, 3.81132e-03, 4.26536e-03, 4.77507e-03, 5.05294e-03, 5.34739e-03], [2.81537e-03, 2.95200e-03, 3.24599e-03, 3.57027e-03, 3.92793e-03, 4.32246e-03, 4.53473e-03, 4.75768e-03, 5.23785e-03, 5.76765e-03, 6.35222e-03, 6.99730e-03, 7.34446e-03, 7.70916e-03, 8.49478e-03, 9.36187e-03, 1.03189e-02, 1.13754e-02, 1.19441e-02, 1.25417e-02], [6.14216e-03, 6.40864e-03, 6.97750e-03, 7.59781e-03, 8.27424e-03, 9.01187e-03, 9.40534e-03, 9.81623e-03, 1.06934e-02, 1.16498e-02, 1.26929e-02, 1.38803e-02, 1.44371e-02, 1.50707e-02, 1.64235e-02, 1.78989e-02, 1.95083e-02, 2.12640e-02, 2.22009e-02, 2.31795e-02], [2.23096e-02, 2.30790e-02, 2.46978e-02, 2.64300e-02, 2.82832e-02, 3.02661e-02, 3.13090e-02, 3.23878e-02, 3.46580e-02, 3.70873e-02, 3.96872e-02, 4.24699e-02, 4.39340e-02, 4.54488e-02, 4.86383e-02, 5.20542e-02, 5.57135e-02, 5.96350e-02, 6.17003e-02, 6.38389e-02], [5.04022e-02, 5.18492e-02, 5.48682e-02, 5.80617e-02, 6.14403e-02, 6.50152e-02, 6.68798e-02, 6.87981e-02, 7.28020e-02, 7.70407e-02, 8.15290e-02, 8.62830e-02, 8.87650e-02, 9.13200e-02, 9.66589e-02, 1.02320e-01, 1.08326e-01, 1.14701e-01, 1.18035e-01, 1.21472e-01], [1.38018e-01, 1.41026e-01, 1.47244e-01, 1.53743e-01, 1.60536e-01, 1.67641e-01, 1.71315e-01, 1.75074e-01, 1.82852e-01, 1.90997e-01, 1.99528e-01, 2.08471e-01, 2.13103e-01, 2.17848e-01, 2.27689e-01, 2.38022e-01, 2.48882e-01, 2.60304e-01, 2.66239e-01, 2.72329e-01], [2.56001e-01, 2.60576e-01, 2.69992e-01, 2.79776e-01, 2.89946e-01, 3.00525e-01, 3.05974e-01, 3.11533e-01, 3.22994e-01, 3.34935e-01, 3.47383e-01, 3.60369e-01, 3.67073e-01, 3.73925e-01, 3.88089e-01, 4.02900e-01, 4.18404e-01, 4.34647e-01, 4.43063e-01, 4.51685e-01], [3.92181e-01, 3.98213e-01, 4.10597e-01, 4.23427e-01, 4.36726e-01, 4.50519e-01, 4.57610e-01, 4.64834e-01, 4.79700e-01, 4.95148e-01, 5.11214e-01, 5.27935e-01, 5.36554e-01, 5.45354e-01, 5.63515e-01, 5.82470e-01, 6.02273e-01, 6.22986e-01, 6.33705e-01, 6.44677e-01], [5.37913e-01, 5.45268e-01, 5.60350e-01, 5.75948e-01, 5.92092e-01, 6.08811e-01, 6.17396e-01, 6.26136e-01, 6.44104e-01, 6.62752e-01, 6.82122e-01, 7.02260e-01, 7.12631e-01, 7.23214e-01, 7.45041e-01, 7.67800e-01, 7.91559e-01, 8.16391e-01, 8.29235e-01, 8.42380e-01], [6.87470e-01, 6.96021e-01, 7.13543e-01, 7.31650e-01, 7.50373e-01, 7.69748e-01, 7.79591e-01, 7.89811e-01, 8.10602e-01, 8.32167e-01, 8.54544e-01, 8.77814e-01, 8.89791e-01, 9.02008e-01, 9.27198e-01, 9.53454e-01, 9.80856e-01, 1.00949e+00, 1.02430e+00, 1.03945e+00], [8.37159e-01, 8.46790e-01, 8.66525e-01, 8.86910e-01, 9.07979e-01, 9.29772e-01, 9.40953e-01, 9.52331e-01, 9.75701e-01, 9.99934e-01, 1.02508e+00, 1.05121e+00, 1.06466e+00, 1.07838e+00, 1.10667e+00, 1.13615e+00, 1.16692e+00, 1.19907e+00, 1.21570e+00, 1.23272e+00], [1.12850e+00, 1.14002e+00, 1.16362e+00, 1.18799e+00, 1.21317e+00, 1.23921e+00, 1.25257e+00, 1.26616e+00, 1.29408e+00, 1.32303e+00, 1.35307e+00, 1.38429e+00, 1.40036e+00, 1.41676e+00, 1.45057e+00, 1.48582e+00, 1.52263e+00, 1.56111e+00, 1.58102e+00, 1.60140e+00], [1.40232e+00, 1.41545e+00, 1.44232e+00, 1.47007e+00, 1.49875e+00, 1.52842e+00, 1.54364e+00, 1.55913e+00, 1.59095e+00, 1.62396e+00, 1.65823e+00, 1.69385e+00, 1.71220e+00, 1.73092e+00, 1.76954e+00, 1.80983e+00, 1.85192e+00, 1.89596e+00, 1.91876e+00, 1.94211e+00], [1.65584e+00, 1.67034e+00, 1.70004e+00, 1.73072e+00, 1.76244e+00, 1.79526e+00, 1.81210e+00, 1.82925e+00, 1.86448e+00, 1.90104e+00, 1.93902e+00, 1.97852e+00, 1.99887e+00, 2.01964e+00, 2.06251e+00, 2.10727e+00, 2.15406e+00, 2.20304e+00, 2.22842e+00, 2.25442e+00], [1.77496e+00, 1.79009e+00, 1.82107e+00, 1.85307e+00, 1.88617e+00, 1.92042e+00, 1.93800e+00, 1.95590e+00, 1.99269e+00, 2.03087e+00, 2.07055e+00, 2.11182e+00, 2.13309e+00, 2.15480e+00, 2.19963e+00, 2.24644e+00, 2.29539e+00, 2.34666e+00, 2.37323e+00, 2.40045e+00], [1.99863e+00, 2.01490e+00, 2.04822e+00, 2.08265e+00, 2.11827e+00, 2.15555e+00, 2.17409e+00, 2.19337e+00, 2.23302e+00, 2.27419e+00, 2.31700e+00, 2.36154e+00, 2.38451e+00, 2.40798e+00, 2.45641e+00, 2.50703e+00, 2.56000e+00, 2.61552e+00, 2.64430e+00, 2.67381e+00], [2.29325e+00, 2.31100e+00, 2.34738e+00, 2.38501e+00, 2.42395e+00, 2.46429e+00, 2.48401e+00, 2.50612e+00, 2.54955e+00, 2.59468e+00, 2.64162e+00, 2.69053e+00, 2.71576e+00, 2.74154e+00, 2.79481e+00, 2.85052e+00, 2.90887e+00, 2.97009e+00, 3.00185e+00, 3.03442e+00], [3.08887e+00, 3.11036e+00, 3.15445e+00, 3.20013e+00, 3.24748e+00, 3.29664e+00, 3.32192e+00, 3.34770e+00, 3.40081e+00, 3.45611e+00, 3.51376e+00, 3.57394e+00, 3.60503e+00, 3.63684e+00, 3.70268e+00, 3.77170e+00, 3.84418e+00, 3.92040e+00, 3.96003e+00, 4.00073e+00], [4.07599e+00, 4.10219e+00, 4.15606e+00, 4.21199e+00, 4.27010e+00, 4.33056e+00, 4.36172e+00, 4.39353e+00, 4.45918e+00, 4.52772e+00, 4.59935e+00, 4.67433e+00, 4.71316e+00, 4.75293e+00, 4.83543e+00, 4.92217e+00, 5.01353e+00, 5.10992e+00, 5.16014e+00, 5.21181e+00], [4.69647e+00, 4.72577e+00, 4.78608e+00, 4.84877e+00, 4.91402e+00, 4.98200e+00, 5.01707e+00, 5.05290e+00, 5.12695e+00, 5.20436e+00, 5.28542e+00, 5.37040e+00, 5.41445e+00, 5.45962e+00, 5.55344e+00, 5.65226e+00, 5.79496e+00, 5.90517e+00, 5.96269e+00, 6.02191e+00], [5.32590e+00, 5.35848e+00, 5.42560e+00, 5.49547e+00, 5.56830e+00, 5.64429e+00, 5.68353e+00, 5.72366e+00, 5.80666e+00, 5.89359e+00, 5.98474e+00, 6.08046e+00, 6.13015e+00, 6.18112e+00, 6.28715e+00, 6.39903e+00, 6.51728e+00, 6.64249e+00, 6.70792e+00, 6.77535e+00]]
# -> Problems for theta < 0.64 corresponding to int(Z) <= 4 ( theta_K_array(Z=5) = 0.64549 and theta_K_array(Z=4) = 0.63640 )
#    NEEDS SMOOTHING IN THIS CASE
	B = maths_library.evaluate_value_from_2D_table(theta_K_array,eta_K_array,B_K_table,theta,eta)
	return max(0.,B);

def eta_K(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : eta = (me/m) E / (Zeff^2 Ry) in ()
#
#   Slater's rule :
	if (int(Z) == 1):
		Zeff = float(Z)
	else :
		Zeff = float(Z) - 0.30
	Zeff2 = Zeff * Zeff
	E     = nrj * 1.e6 / Ry
	a     = me / m
	eta   = a * E / Zeff2
	return eta;

def C_K(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C     = K-shell correction contribution to the stopping number in ()
#
	theta    = theta_K(Z)
	S        = S_K(Z)
	T        = T_K(Z)
	U        = U_K(Z)
	V        = V_K(Z)
	eta      = eta_K(Z,m,nrj)
	ieta     = 1. / eta
	ieta_min = 0.08
	ieta_max = 0.12
	if (ieta < ieta_min) :
#   according to M. C. Walske, Phys. Rev. Vol. 88 No. 6 P. 1283 (1952):
#   C_K(Z,eta_K) = (U_K / eta_K) + (V_K / eta_K^2)
		ieta2 = ieta * ieta
		C = (U*ieta) + (V*ieta2)
	elif (ieta >= ieta_min) and (ieta < ieta_max):
#   linear interpolation between the two expressions
		ieta_min2 = ieta_min * ieta_min
		Cmin      = (U*ieta_min) + (V*ieta_min2)
		Bmax      = B_K(theta,1./ieta_max)
		Cmax      = ( S * math.log(1./ieta_max) ) + T - Bmax
		a         = (Cmax - Cmin) / (ieta_max - ieta_min)
		b         = Cmax - ( a * ieta_max )
		C         = (a * ieta) + b
	else :
#   according to M. S. Livingston and H. A. Bethe, Revs. Modern Phys. 9, 263 (1937) :
#   C_K(Z,eta_K) = S(Z) log(eta_K) + T(Z) - B_K(Z,eta_K)
		B = B_K(theta,eta)
		C = ( S * math.log(eta) ) + T - B
	return C;

def theta_L(Z):
# input  : Z     = Atomic number in ()
# output : theta = 1 - (I_R - I_0 )/I_NR where I_NR = theoretical non-relativistic ionization potential of L-shell electron(s)
#                                              I_R  = theoretical relativistic ionization potential of L-shell electron(s)
#                                              I_0  = experimental energy difference from the L-shell to the lowest-lying unoccupied state
#   M. C. Walske, Phys. Rev., Vol. 101, No. 3, p. 940 (1956) and
#   V. H. Hönl, Z. Physik 84, 1 (1933) :
	theta_L_array = [0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.     , 0.04235683, 0.11040844, 0.1653322, 0.21049141, 0.24821736, 0.28014611, 0.30769309, 0.3314791, 0.35232331, 0.37174392, 0.38267   , 0.39281   , 0.40468  , 0.41483  , 0.42273   , 0.43227  , 0.44185   , 0.45037   , 0.45882   , 0.46626   , 0.47308   , 0.48034   , 0.48766   , 0.49452   , 0.50083   , 0.50702   , 0.51236   , 0.51826   , 0.52282  , 0.52676   , 0.53109   , 0.53521   , 0.53879   , 0.54217   , 0.54588   , 0.54916   , 0.55197   , 0.55510   , 0.55787   , 0.56070   , 0.56379   , 0.56634   , 0.56884  , 0.57165  , 0.57446   , 0.57727   , 0.57966   , 0.58155   , 0.58418   , 0.58701   , 0.58886   , 0.59127   , 0.59393   , 0.59604   , 0.59785   , 0.59987  , 0.60109   , 0.60287   , 0.60503   , 0.60671   , 0.60840   , 0.61004   , 0.61128   , 0.61399   , 0.61513   , 0.61681   , 0.61793   , 0.62017   , 0.62187   , 0.62409   , 0.62475]	
	z             = int(Z) - 1
	theta         = theta_L_array[z]
	return theta;

def S_L(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the L-shell contribution to the stopping number
#              B_L(Z,eta_K) = S(Z) log(eta_L) + T(Z) - C_L(Z,eta_L)
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	S_L_array = [32.67685, 29.25826, 26.48722, 24.19565, 22.26903, 20.62660, 19.20980, 17.97512, 16.88957, 15.92767, 14.94847, 14.02541, 13.25067, 12.54122, 11.90843, 11.41276, 10.89866, 10.45111, 10.06763, 9.70581, 9.43481, 9.18505, 8.96004, 8.70945, 8.51198, 8.36456, 8.20090, 8.03932, 7.90623, 7.78151, 7.67927, 7.58686, 7.48897, 7.40022, 7.31704, 7.24155, 7.17411, 7.11593, 7.05165, 7.00499, 6.96628, 6.92374, 6.88326, 6.84808, 6.81645, 6.78269, 6.75284, 6.72810, 6.70093, 6.67688, 6.65274, 6.62779, 6.60720, 6.58701, 6.56432, 6.54163, 6.51894, 6.49964, 6.48549, 6.46615, 6.44534, 6.43173, 6.41400, 6.39444, 6.37892, 6.36561, 6.35075, 6.34247, 6.33051, 6.31599, 6.30470, 6.29335, 6.28233, 6.27399, 6.25578, 6.24812, 6.23683, 6.22931, 6.21435, 6.20389, 6.19022, 6.18616]
	z         = int(Z) - 1
	S         = S_L_array[z]
	return S; 

def T_L(Z):
# input  : Z = Atomic number in ()
# output : T = T coefficient of asymptotic stopping number formula for the L-shell contribution to the stopping number
#              B_K(Z,eta_L) = S(Z) log(eta_L) + T(Z) - C_K(Z,eta_L)
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	T_L_array = [51.26351, 48.50790, 46.15369, 44.11202, 42.31934, 40.72884, 39.30517, 38.02105, 36.85510, 35.79023, 34.61551, 33.53561, 32.57885, 31.67825, 30.84430, 30.15874, 29.43174, 28.77285, 28.19169, 27.62092, 27.17845, 26.76576, 26.37978, 25.94222, 25.58682, 25.31747, 25.00880, 24.70216, 24.44219, 24.19308, 23.98330, 23.79264, 23.59032, 23.39923, 23.22015, 23.05680, 22.90542, 22.77483, 22.63055, 22.52321, 22.43268, 22.33320, 22.23855, 22.15629, 22.08086, 21.99943, 21.92743, 21.86699, 21.80026, 21.74121, 21.68149, 21.61837, 21.56627, 21.51519, 21.45779, 21.40038, 21.34297, 21.29414, 21.25723, 21.20640, 21.15169, 21.11593, 21.06935, 21.01793, 20.97714, 20.94215, 20.90311, 20.88062, 20.84800, 20.80842, 20.77763, 20.74666, 20.71661, 20.69389, 20.64423, 20.62334, 20.59255, 20.57203, 20.53114, 20.50154, 20.46289, 20.45140]
	z         = int(Z) - 1
	T         = T_L_array[z]
	return T; 

def U_L(Z):
# input  : Z = Atomic number in ()
# output : U = U coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              C_L(Z,eta_K) = (U_L / eta_K)
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	U_L_array = [-3.06925, -2.73850, -2.40775, -2.07700, -1.74625, -1.41549, -1.08475, -0.754, -0.42325, -0.09250, 0.23348, 0.50138, 0.7151, 0.90548, 1.06868, 1.18978, 1.31211, 1.41331, 1.49667, 1.57101, 1.62373, 1.67138, 1.71159, 1.75495, 1.78718, 1.81052, 1.83470, 1.85824, 1.87629, 1.89226, 1.90442, 1.91523, 1.92661, 1.93562, 1.94405, 1.95158, 1.95739, 1.96241, 1.96796, 1.97157, 1.97433, 1.97736, 1.98024, 1.98275, 1.98477, 1.98677, 1.98854, 1.98988, 1.99129, 1.99254, 1.99373, 1.99475, 1.99559, 1.99641, 1.99734, 1.99827, 1.99919, 1.99998, 2.00039, 2.00089, 2.00143, 2.00178, 2.00224, 2.00274, 2.00314, 2.00349, 2.00387, 2.00397, 2.00410, 2.00425, 2.00436, 2.00448, 2.00460, 2.00468, 2.00487, 2.00495, 2.00507, 2.00515, 2.00529, 2.00526, 2.00521, 2.00520]
#   The extrapolation gives U_L < 0 for int(Z) < 11!!!
	if (int(Z) < 11):
		z = 9
	else:
		z = int(Z) - 1
	U = U_L_array[z]
	return U; 

def B_LI(theta,eta):
# input  : theta = 1 - (I_R - I_0 )/I_NR where I_NR = theoretical non-relativistic ionization potential of L-shell electron(s)
#                                              I_R  = theoretical relativistic ionization potential of L-shell electron(s)
#                                              I_0  = experimental energy difference from the L-shell to the lowest-lying unoccupied state
#          eta   = me v^2 / (2 Zeff^2 Ry) where v is the projectile velocity and Zeff the effective charge seen by a LI-shell electron
# output : B     = LI-shell contribution to the stopping number
#
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 : B_L_table[i][j] = B_L(eta_L[i],theta_L[j]) with :
	theta_LI_array = [0.66, 0.65, 0.64, 0.62, 0.6, 0.58, 0.56, 0.55, 0.5, 0.52, 0.5, 0.48, 0.46, 0.45, 0.44, 0.42, 0.4, 0.38, 0.36, 0.35, 0.34, 0.32, 0.3, 0.28, 0.26, 0.24]
	eta_LI_array   = [0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.5, 1.7, 2.0, 3.0, 5.0, 7.0]
	B_LI_table     = [[2.41110e-04, 2.56120e-04, 2.72020e-04, 3.06580e-04, 3.45110e-04, 3.87950e-04, 4.35420e-04, 4.61000e-04, 4.87860e-04, 5.45610e-04, 6.09050e-04, 6.78630e-04, 7.54940e-04, 7.95870e-04, 8.38830e-04, 9.31600e-04, 1.03520e-03, 1.15290e-03, 1.28950e-03, 1.36700e-03, 1.45240e-03, 1.65240e-03, 1.90780e-03, 2.24140e-03, 2.68890e-03, 3.30060e-03], [6.39470e-04, 6.70580e-04, 7.02950e-04, 7.71670e-04, 8.45920e-04, 9.26050e-04, 1.01250e-03, 1.05830e-03, 1.10580e-03, 1.20680e-03, 1.31700e-03, 1.43770e-03, 1.57190e-03, 1.64510e-03, 1.72310e-03, 1.89690e-03, 2.10090e-03, 2.34590e-03, 2.64670e-03, 2.82420e-03, 3.02380e-03, 3.50450e-03, 4.12600e-03, 4.93760e-03, 6.00500e-03, 7.41520e-03], [1.54690e-03, 1.60360e-03, 1.66220e-03, 1.78560e-03, 1.91810e-03, 2.16150e-03, 2.31780e-03, 2.40190e-03, 2.49040e-03, 2.68320e-03, 2.90170e-03, 3.15340e-03, 3.44790e-03, 3.61490e-03, 3.79760e-03, 4.21870e-03, 4.73200e-03, 5.36360e-03, 6.14720e-03, 6.60860e-03, 7.12460e-03, 8.34910e-03, 9.88710e-03, 1.18220e-02, 1.42610e-02, 1.73350e-02], [3.72210e-03, 3.84040e-03, 3.96500e-03, 4.23540e-03, 4.53960e-03, 4.88530e-03, 5.28200e-03, 5.50310e-03, 5.74140e-03, 6.27750e-03, 6.90770e-03, 7.65250e-03, 8.53700e-02, 9.04070e-03, 9.59100e-03, 1.08500e-02, 1.23580e-02, 1.41650e-02, 1.63316e-02, 1.75720e-02, 1.89320e-02, 2.20530e-02, 2.58030e-02, 3.03080e-02, 3.57280e-02, 4.22580e-02], [6.94490e-03, 7.19100e-03, 7.45350e-03, 8.03360e-03, 8.69840e-03, 9.46380e-03, 1.03480e-02, 1.08410e-02, 1.13720e-02, 1.25610e-02, 1.39430e-02, 1.55530e-02, 1.73270e-02, 1.84780e-02, 1.96120e-02, 2.21600e-02, 2.51300e-02, 2.85940e-02, 3.26340e-02, 3.49000e-02, 3.73480e-02, 4.28500e-02, 4.92780e-02, 5.67980e-02, 6.56100e-02, 7.59630e-02], [1.74110e-02, 1.81800e-02, 1.89970e-02, 2.07840e-02, 2.27970e-02, 2.50660e-02, 2.76220e-02, 2.90200e-02, 3.05030e-02, 3.37500e-02, 3.74700e-02, 4.15280e-02, 4.61700e-02, 4.87080e-02, 5.14010e-02, 5.72970e-02, 6.39430e-02, 7.14410e-02, 7.99070e-02, 8.45440e-02, 8.94860e-02, 1.00320e-01, 1.12600e-01, 1.26560e-01, 1.42480e-01, 1.60710e-01], [3.84740e-02, 4.00560e-02, 4.17180e-02, 4.53000e-02, 4.92540e-02, 5.36190e-02, 5.84360e-02, 6.10280e-02, 6.37520e-02, 6.96190e-02, 7.60980e-02, 8.32490e-02, 9.11500e-02, 9.54060e-02, 9.98810e-02, 1.09540e-01, 1.20230e-01, 1.32080e-01, 1.45230e-01, 1.52370e-01, 1.59850e-01, 1.76140e-01, 1.94340e-01, 2.14730e-01, 2.37660e-01, 2.63570e-01], [6.73710e-02, 6.99280e-02, 7.25960e-02, 7.82820e-02, 8.44700e-02, 9.12060e-02, 9.85380e-02, 1.02440e-01, 1.06520e-01, 1.15220e-01, 1.24700e-01, 1.35040e-01, 1.46320e-01, 1.52340e-01, 1.58640e-01, 1.72110e-01, 1.86860e-01, 2.03040e-01, 2.20820e-01, 2.30360e-01, 2.40380e-01, 2.61990e-01, 2.85900e-01, 3.12480e-01, 3.42120e-01, 3.75360e-01], [1.04180e-01, 1.07780e-01, 1.11520e-01, 1.19430e-01, 1.27960e-01, 1.37150e-01, 1.47060e-01, 1.52310e-01, 1.57760e-01, 1.69310e-01, 1.81790e-01, 1.95300e-01, 2.09910e-01, 2.17670e-01, 2.25760e-01, 2.42950e-01, 2.61650e-01, 2.82010e-01, 3.04230e-01, 3.16110e-01, 3.28540e-01, 3.55220e-01, 3.84590e-01, 4.17040e-01, 4.53060e-01, 4.93260e-01], [1.96470e-01, 2.02170e-01, 2.08050e-01, 2.20380e-01, 2.33510e-01, 2.47510e-01, 2.62440e-01, 2.70270e-01, 2.78370e-01, 2.95400e-01, 3.13610e-01, 3.33120e-01, 3.54030e-01, 3.65060e-01, 3.76500e-01, 4.00670e-01, 4.26730e-01, 4.54880e-01, 4.85360e-01, 5.01560e-01, 5.18460e-01, 5.54530e-01, 5.93970e-01, 6.37280e-01, 6.85070e-01, 7.38100e-01], [3.05940e-01, 3.13610e-01, 3.21500e-01, 3.37960e-01, 3.55370e-01, 3.73810e-01, 3.93360e-01, 4.03570e-01, 4.14100e-01, 4.36130e-01, 4.59560e-01, 4.85200e-01, 5.11150e-01, 5.25140e-01, 5.39610e-01, 5.70080e-01, 6.02770e-01, 6.37930e-01, 6.75860e-01, 6.95960e-01, 7.16880e-01, 7.61410e-01, 8.09920e-01, 8.63010e-01, 9.21420e-01, 9.86040e-01], [6.14110e-01, 6.25970e-01, 6.38110e-01, 6.63300e-01, 6.89740e-01, 7.17530e-01, 7.46780e-01, 7.61990e-01, 7.77610e-01, 8.10140e-01, 8.44530e-01, 8.80930e-01, 9.19540e-01, 9.39730e-01, 9.60560e-01, 1.00430e+00, 1.05090e+00, 1.10080e+00, 1.15440e+00, 1.18280e+00, 1.21220e+00, 1.27460e+00, 1.34240e+00, 1.41630e+00, 1.49740e+00, 1.58680e+00], [9.31000e-01, 9.56140e-01, 9.71620e-01, 1.00370e+00, 1.03720e+00, 1.07230e+00, 1.10920e+00, 1.12840e+00, 1.14800e+00, 1.18880e+00, 1.23190e+00, 1.27740e+00, 1.32550e+00, 1.35060e+00, 1.37650e+00, 1.43080e+00, 1.48860e+00, 1.55040e+00, 1.61670e+00, 1.65170e+00, 1.68800e+00, 1.76500e+00, 1.84840e+00, 1.93940e+00, 2.03900e+00, 2.14890e+00], [1.51720e+00, 1.53720e+00, 1.55760e+00, 1.59980e+00, 1.64380e+00, 1.68990e+00, 1.73820e+00, 1.76320e+00, 1.78890e+00, 1.84220e+00, 1.89830e+00, 1.95750e+00, 2.02010e+00, 2.05280e+00, 2.08640e+00, 2.15690e+00, 2.23190e+00, 2.31200e+00, 2.39790e+00, 2.44320e+00, 2.49020e+00, 2.58990e+00, 2.69800e+00, 2.81590e+00, 2.94510e+00, 3.08760e+00], [2.01730e+00, 2.04080e+00, 2.06470e+00, 2.11420e+00, 2.16590e+00, 2.21990e+00, 2.27650e+00, 2.30590e+00, 2.33600e+00, 2.39840e+00, 2.46420e+00, 2.53360e+00, 2.60700e+00, 2.64520e+00, 2.68470e+00, 2.76740e+00, 2.85540e+00, 2.94940e+00, 3.05020e+00, 3.10340e+00, 3.15860e+00, 3.27580e+00, 3.40300e+00, 3.54160e+00, 3.69380e+00, 3.86200e+00], [2.39320e+00, 2.41930e+00, 2.44600e+00, 2.50110e+00, 2.55870e+00, 2.61900e+00, 2.68210e+00, 2.71480e+00, 2.74840e+00, 2.81810e+00, 2.89150e+00, 2.96900e+00, 3.05090e+00, 3.09370e+00, 3.13780e+00, 3.23010e+00, 3.32850e+00, 3.43370e+00, 3.54660e+00, 3.60620e+00, 3.66810e+00, 3.79940e+00, 3.94210e+00, 4.09780e+00, 4.26880e+00, 4.45800e+00], [2.70910e+00, 2.73740e+00, 2.76630e+00, 2.82600e+00, 2.88840e+00, 2.95380e+00, 3.02220e+00, 3.05770e+00, 3.09410e+00, 3.16980e+00, 3.24940e+00, 3.33360e+00, 3.42260e+00, 3.46910e+00, 3.51710e+00, 3.61750e+00, 3.72460e+00, 3.83910e+00, 3.96200e+00, 4.02700e+00, 4.09450e+00, 4.23780e+00, 4.39350e+00, 4.56360e+00, 4.75060e+00, 4.95760e+00], [2.97420e+00, 3.00440e+00, 3.03520e+00, 3.09880e+00, 3.16520e+00, 3.23490e+00, 3.30790e+00, 3.34570e+00, 3.38450e+00, 3.46520e+00, 3.55020e+00, 3.64000e+00, 3.73510e+00, 3.78480e+00, 3.83600e+00, 3.94330e+00, 4.05780e+00, 4.18040e+00, 4.30200e+00, 4.37150e+00, 4.44380e+00, 4.59740e+00, 4.76440e+00, 4.94700e+00, 5.14780e+00, 5.37030e+00], [3.22220e+00, 3.25390e+00, 3.28630e+00, 3.35320e+00, 3.42320e+00, 3.49650e+00, 3.57340e+00, 3.61330e+00, 3.65420e+00, 3.73920e+00, 3.82890e+00, 3.92360e+00, 4.02390e+00, 4.07640e+00, 4.13040e+00, 4.24380e+00, 4.36480e+00, 4.49440e+00, 4.63360e+00, 4.70720e+00, 4.78370e+00, 4.94630e+00, 5.12330e+00, 5.31690e+00, 5.53000e+00, 5.76610e+00], [3.66900e+00, 3.70330e+00, 3.73840e+00, 3.81080e+00, 3.88660e+00, 3.96610e+00, 4.04950e+00, 4.09280e+00, 4.13710e+00, 4.22950e+00, 4.32690e+00, 4.42990e+00, 4.53910e+00, 4.59620e+00, 4.65510e+00, 4.77860e+00, 4.91060e+00, 5.05200e+00, 5.20410e+00, 5.28450e+00, 5.36820e+00, 5.54620e+00, 5.74000e+00, 5.95230e+00, 6.18630e+00, 6.44580e+00], [3.98190e+00, 4.01830e+00, 4.05550e+00, 4.13240e+00, 4.21300e+00, 4.29740e+00, 4.38610e+00, 4.43210e+00, 4.47940e+00, 4.57770e+00, 4.68140e+00, 4.79120e+00, 4.90760e+00, 4.96850e+00, 5.03140e+00, 5.16330e+00, 5.30430e+00, 5.45550e+00, 5.61820e+00, 5.70440e+00, 5.79400e+00, 5.98480e+00, 6.19270e+00, 6.42060e+00, 6.67190e+00, 6.95100e+00], [4.27450e+00, 4.31270e+00, 4.35170e+00, 4.43240e+00, 4.51700e+00, 4.60560e+00, 4.69880e+00, 4.74710e+00, 4.79680e+00, 4.90010e+00, 5.00920e+00, 5.12470e+00, 5.24730e+00, 5.31140e+00, 5.37760e+00, 5.51660e+00, 5.66530e+00, 5.82490e+00, 5.99670e+00, 6.08760e+00, 6.18230e+00, 6.38390e+00, 6.60380e+00, 6.84510e+00, 7.11130e+00, 7.40710e+00], [4.40470e+00, 4.44360e+00, 4.48340e+00, 4.56580e+00, 4.65210e+00, 4.74260e+00, 4.83780e+00, 4.88720e+00, 4.93790e+00, 5.04340e+00, 5.15500e+00, 5.27310e+00, 5.39840e+00, 5.46400e+00, 5.53170e+00, 5.67390e+00, 5.82600e+00, 5.98930e+00, 6.16520e+00, 6.25830e+00, 6.35530e+00, 6.56180e+00, 6.78710e+00, 7.03440e+00, 7.30730e+00, 7.61070e+00], [4.63830e+00, 4.67870e+00, 4.72000e+00, 4.80540e+00, 4.89490e+00, 4.98880e+00, 5.08760e+00, 5.13880e+00, 5.19150e+00, 5.30110e+00, 5.41700e+00, 5.53980e+00, 5.67010e+00, 5.73730e+00, 5.80880e+00, 5.95680e+00, 6.11520e+00, 6.28530e+00, 6.46860e+00, 6.56570e+00, 6.66680e+00, 6.88230e+00, 7.11750e+00, 7.37570e+00, 7.66090e+00, 7.97820e+00], [4.93690e+00, 4.97910e+00, 5.02230e+00, 5.11160e+00, 5.20530e+00, 5.30360e+00, 5.40700e+00, 5.46070e+00, 5.51590e+00, 5.63080e+00, 5.75230e+00, 5.88110e+00, 6.01740e+00, 6.08960e+00, 6.16360e+00, 6.31920e+00, 6.48570e+00, 6.66470e+00, 6.85770e+00, 6.96000e+00, 7.06660e+00, 7.29370e+00, 7.54180e+00, 7.81430e+00, 8.11560e+00, 8.45100e+00], [5.65140e+00, 5.69810e+00, 5.74590e+00, 5.84500e+00, 5.94890e+00, 6.05810e+00, 6.17300e+00, 6.23280e+00, 6.29430e+00, 6.42240e+00, 6.55800e+00, 6.70190e+00, 6.85490e+00, 6.93510e+00, 7.01800e+00, 7.19250e+00, 7.37950e+00, 7.58080e+00, 7.79810e+00, 7.91340e+00, 8.03360e+00, 8.29010e+00, 8.57080e+00, 8.87960e+00, 9.22140e+00, 9.60270e+00], [6.46650e+00, 6.51890e+00, 6.57240e+00, 6.68350e+00, 6.80030e+00, 6.92310e+00, 7.05250e+00, 7.11990e+00, 7.18920e+00, 7.33380e+00, 7.48720e+00, 7.65000e+00, 7.82350e+00, 7.91460e+00, 8.00870e+00, 8.20710e+00, 8.42000e+00, 8.64960e+00, 8.89780e+00, 9.02970e+00, 9.16730e+00, 9.46120e+00, 9.78340e+00, 1.01384e+01, 1.05323e+01, 1.09722e+01], [6.86340e+00, 6.91940e+00, 6.97670e+00, 7.09570e+00, 7.22080e+00, 7.35260e+00, 7.49150e+00, 7.56390e+00, 7.63840e+00, 7.79380e+00, 7.95880e+00, 8.13420e+00, 8.32110e+00, 8.41930e+00, 8.52090e+00, 8.73500e+00, 8.96510e+00, 9.21330e+00, 9.48190e+00, 9.62480e+00, 9.77390e+00, 1.00926e+01, 1.04423e+01, 1.08282e+01, 1.12565e+01, 1.17356e+01]]
# -> Problems for theta < 0.24 corresponding to int(Z) <= 10 ( theta_L_array(Z=11) = 0.24553 and theta_L_array(Z=10) = 0.23230 )
#    NEEDS SMOOTHING IN THIS CASE
	B = maths_library.evaluate_value_from_2D_table(theta_LI_array,eta_LI_array,B_LI_table,theta,eta)
	return max(0.,B);

def B_LII(theta,eta):
# input  : theta = 1 - (I_R - I_0 )/I_NR where I_NR = theoretical non-relativistic ionization potential of L-shell electron(s)
#                                              I_R  = theoretical relativistic ionization potential of L-shell electron(s)
#                                              I_0  = experimental energy difference from the L-shell to the lowest-lying unoccupied state
#          eta   = me v^2 / (2 Zeff^2 Ry) where v is the projectile velocity and Zeff the effective charge seen by a LII-shell electron
# output : B     = LII-shell contribution to the stopping number
#
#   according to G. S. Khandelwal, Nucl. Phys. A116 (1966) p. 97-111 :
	theta_LII_array = [0.66, 0.65, 0.64, 0.62, 0.6, 0.58, 0.56, 0.55, 0.54, 0.52, 0.5, 0.48, 0.46, 0.45, 0.44, 0.42, 0.4, 0.38, 0.36, 0.35, 0.34, 0.32, 0.3, 0.28, 0.26, 0.24]
	eta_LII_array   = [0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.5, 1.7, 2.0, 3.0, 5.0, 7.0]
	B_LII_table     = [[3.63240e-05, 4.06090e-05, 4.54300e-05, 5.69690e-05, 7.16250e-05, 9.02790e-05, 1.14070e-04, 1.28340e-04, 1.44470e-04, 1.83390e-04, 2.33300e-04, 2.97380e-04, 3.79770e-04, 4.29450e-04, 4.85820e-04, 6.22440e-04, 7.98580e-04, 1.02580e-03, 1.31900e-03, 1.49610e-03, 1.69740e-03, 2.18580e-03, 2.81630e-03, 3.63020e-03, 4.68140e-03, 6.03950e-03], [1.81100e-04, 2.00010e-04, 2.20990e-04, 2.70060e-04, 3.30490e-04, 4.04980e-04, 4.96880e-04, 5.50610e-04, 6.10320e-04, 7.50420e-04, 9.23550e-04, 1.13750e-03, 1.40210e-03, 1.55700e-03, 1.72920e-03, 2.13350e-03, 2.63350e-03, 3.25150e-03, 4.01580e-03, 4.46230e-03, 4.95920e-03, 6.12570e-03, 7.56750e-03, 9.35020e-03, 1.15560e-02, 1.42900e-02], [8.65240e-04, 9.42230e-04, 1.02620e-03, 1.21780e-03, 1.44590e-03, 1.71740e-03, 2.04050e-03, 2.22450e-03, 2.42520e-03, 2.88290e-03, 3.42750e-03, 4.07580e-03, 4.84570e-03, 5.28390e-03, 5.76170e-03, 6.85040e-03, 8.14420e-03, 9.68160e-03, 1.15090e-02, 1.25480e-02, 1.36810e-02, 1.62630e-02, 1.93360e-02, 2.29990e-02, 2.73700e-02, 3.26030e-02], [4.22930e-03, 4.53140e-03, 4.85510e-03, 5.57310e-03, 6.39680e-03, 7.34140e-03, 8.42420e-03, 9.02360e-03, 9.66520e-03, 1.10870e-02, 1.27160e-02, 1.45810e-02, 1.67170e-02, 1.78980e-02, 1.91630e-02, 2.19640e-02, 2.51730e-02, 2.88510e-02, 3.30700e-02, 3.54080e-02, 3.79140e-02, 4.34830e-02, 4.98980e-02, 5.73040e-02, 6.58840e-02, 7.58610e-02], [1.14850e-02, 1.21720e-02, 1.29000e-02, 1.44860e-02, 1.62640e-02, 1.82560e-02, 2.04870e-02, 2.17020e-02, 2.29890e-02, 2.57860e-02, 2.89220e-02, 3.24350e-02, 3.63710e-02, 3.85140e-02, 4.07840e-02, 4.57330e-02, 5.12880e-02, 5.75310e-02, 6.45550e-02, 6.83940e-02, 7.24720e-02, 8.14130e-02, 9.15390e-02, 1.03040e-01, 1.16170e-01, 1.31210e-01], [3.94710e-02, 4.12700e-02, 4.31490e-02, 4.71630e-02, 5.15430e-02, 5.64230e-02, 6.15400e-02, 6.43260e-02, 6.72370e-02, 7.34610e-02, 8.02640e-02, 8.77050e-02, 9.58520e-02, 1.00210e-01, 1.04780e-01, 1.14580e-01, 1.25350e-01, 1.37210e-01, 1.50300e-01, 1.51010e-01, 1.58440e-01, 1.74510e-01, 1.92440e-01, 2.12440e-01, 2.34960e-01, 2.60440e-01], [8.44540e-02, 8.75990e-02, 9.08600e-02, 9.77470e-02, 1.05160e-01, 1.13130e-01, 1.21710e-01, 1.26250e-01, 1.30960e-01, 1.40940e-01, 1.51720e-01, 1.63360e-01, 1.75960e-01, 1.82650e-01, 1.89620e-01, 2.04450e-01, 2.20580e-01, 2.38180e-01, 2.57430e-01, 2.67740e-01, 2.78550e-01, 3.01800e-01, 3.27510e-01, 3.56080e-01, 3.88030e-01, 4.24010e-01], [1.43390e-01, 1.47950e-01, 1.52660e-01, 1.62530e-01, 1.73060e-01, 1.84300e-01, 1.96300e-01, 2.02610e-01, 2.09240e-01, 2.22890e-01, 2.37620e-01, 2.53440e-01, 2.70450e-01, 2.79440e-01, 2.88770e-01, 3.08550e-01, 3.29950e-01, 3.53180e-01, 3.78460e-01, 3.91950e-01, 4.06070e-01, 4.36350e-01, 4.69730e-01, 5.06720e-01, 5.47980e-01, 5.94360e-01], [2.13040e-01, 2.18990e-01, 2.25120e-01, 2.37940e-01, 2.51530e-01, 2.65960e-01, 2.81300e-01, 2.89340e-01, 2.97630e-01, 3.15030e-01, 3.33610e-01, 3.53460e-01, 3.74730e-01, 3.85940e-01, 3.97560e-01, 4.22120e-01, 4.48610e-01, 4.77270e-01, 5.08390e-01, 5.24970e-01, 5.42300e-01, 5.79430e-01, 6.20280e-01, 6.65490e-01, 7.15890e-01, 7.72520e-01], [3.73820e-01, 3.82410e-01, 3.91220e-01, 4.09550e-01, 4.28880e-01, 4.49280e-01, 4.70860e-01, 4.82120e-01, 4.93710e-01, 5.17930e-01, 5.43680e-01, 5.71090e-01, 6.00320e-01, 6.15690e-01, 6.31590e-01, 6.65120e-01, 7.01190e-01, 7.40120e-01, 7.82300e-01, 8.04740e-01, 8.28180e-01, 8.78360e-01, 9.33550e-01, 9.94620e-01, 1.06270e+00, 1.13940e+00], [5.50560e-01, 5.61510e-01, 5.72730e-01, 5.96010e-01, 6.20490e-01, 6.46270e-01, 6.73460e-01, 6.87620e-01, 7.02180e-01, 7.32580e-01, 7.64810e-01, 7.99070e-01, 8.35560e-01, 8.54720e-01, 8.74540e-01, 9.16300e-01, 9.61190e-01, 1.00960e+00, 1.06210e+00, 1.09000e+00, 1.11920e+00, 1.18160e+00, 1.25030e+00, 1.32650e+00, 1.41160e+00, 1.50760e+00], [1.00660e+00, 1.02240e+00, 1.03860e+00, 1.07210e+00, 1.10730e+00, 1.14430e+00, 1.18320e+00, 1.20350e+00, 1.22430e+00, 1.26770e+00, 1.31370e+00, 1.36260e+00, 1.41470e+00, 1.44210e+00, 1.47030e+00, 1.53000e+00, 1.59420e+00, 1.66360e+00, 1.73890e+00, 1.77900e+00, 1.82100e+00, 1.91120e+00, 2.01080e+00, 2.12170e+00, 2.24620e+00, 2.38760e+00], [1.43760e+00, 1.45720e+00, 1.47730e+00, 1.51880e+00, 1.56240e+00, 1.60830e+00, 1.65660e+00, 1.68170e+00, 1.70760e+00, 1.76150e+00, 1.81880e+00, 1.87970e+00, 1.94460e+00, 1.97880e+00, 2.01420e+00, 2.08890e+00, 2.16940e+00, 2.25670e+00, 2.35160e+00, 2.40240e+00, 2.45560e+00, 2.57010e+00, 2.69710e+00, 2.83910e+00, 2.99940e+00, 3.18220e+00], [2.17120e+00, 2.19640e+00, 2.22230e+00, 2.27580e+00, 2.33220e+00, 2.39150e+00, 2.45420e+00, 2.48690e+00, 2.52050e+00, 2.59090e+00, 2.66580e+00, 2.74570e+00, 2.83120e+00, 2.87620e+00, 2.92310e+00, 3.02220e+00, 3.12950e+00, 3.24630e+00, 3.37410e+00, 3.44270e+00, 3.51480e+00, 3.67060e+00, 3.84450e+00, 4.04040e+00, 4.26310e+00, 4.51930e+00], [2.75000e+00, 2.77930e+00, 2.80940e+00, 2.87190e+00, 2.93770e+00, 3.00720e+00, 3.08070e+00, 3.11920e+00, 3.15870e+00, 3.24170e+00, 3.33020e+00, 3.42490e+00, 3.52660e+00, 3.58030e+00, 3.63610e+00, 3.75460e+00, 3.88350e+00, 4.02420e+00, 4.17880e+00, 4.26200e+00, 4.34960e+00, 4.53980e+00, 4.75300e+00, 4.99440e+00, 5.27030e+00, 5.58950e+00], [3.20330e+00, 3.23590e+00, 3.26930e+00, 3.33890e+00, 3.41220e+00, 3.48980e+00, 3.57210e+00, 3.61510e+00, 3.65950e+00, 3.75270e+00, 3.85230e+00, 3.95910e+00, 4.07410e+00, 4.13500e+00, 4.19830e+00, 4.33300e+00, 4.47990e+00, 4.64080e+00, 4.81800e+00, 4.91370e+00, 5.01460e+00, 5.23410e+00, 5.48110e+00, 5.76180e+00, 6.08400e+00, 6.45830e+00], [3.60380e+00, 3.63910e+00, 3.67530e+00, 3.75060e+00, 3.83030e+00, 3.91460e+00, 4.00420e+00, 4.05110e+00, 4.09950e+00, 4.20130e+00, 4.31030e+00, 4.42740e+00, 4.55370e+00, 4.62060e+00, 4.69040e+00, 4.83900e+00, 5.00130e+00, 5.17960e+00, 5.37650e+00, 5.48300e+00, 5.59540e+00, 5.84060e+00, 6.11730e+00, 6.43260e+00, 6.79580e+00, 7.21910e+00], [3.91060e+00, 3.94820e+00, 3.98670e+00, 4.06700e+00, 4.15200e+00, 4.24210e+00, 4.33800e+00, 4.38820e+00, 4.44010e+00, 4.54930e+00, 4.66640e+00, 4.79250e+00, 4.92860e+00, 5.00090e+00, 5.07620e+00, 5.23700e+00, 5.41290e+00, 5.60660e+00, 5.82080e+00, 5.93690e+00, 6.05960e+00, 6.32750e+00, 6.63060e+00, 6.97690e+00, 7.37670e+00, 7.84400e+00], [4.17900e+00, 4.21850e+00, 4.25900e+00, 4.34370e+00, 4.43330e+00, 4.52850e+00, 4.62980e+00, 4.68300e+00, 4.73800e+00, 4.85370e+00, 4.97800e+00, 5.11190e+00, 5.25680e+00, 5.33380e+00, 5.41410e+00, 5.58570e+00, 5.77380e+00, 5.98110e+00, 6.21090e+00, 6.33550e+00, 6.46740e+00, 6.75580e+00, 7.08270e+00, 7.45700e+00, 7.89000e+00, 8.39720e+00], [4.63440e+00, 4.67720e+00, 4.72120e+00, 4.81310e+00, 4.91060e+00, 5.01440e+00, 5.12500e+00, 5.18310e+00, 5.24320e+00, 5.37010e+00, 5.50660e+00, 5.65400e+00, 5.81380e+00, 5.89890e+00, 5.98780e+00, 6.17800e+00, 6.38700e+00, 6.61790e+00, 6.87470e+00, 7.01420e+00, 7.16210e+00, 7.48610e+00, 7.85460e+00, 8.27780e+00, 8.76900e+00, 9.34640e+00], [4.97870e+00, 5.02420e+00, 5.07110e+00, 5.16890e+00, 5.27290e+00, 5.38370e+00, 5.50500e+00, 5.56420e+00, 5.62870e+00, 5.76480e+00, 5.91140e+00, 6.07010e+00, 6.24240e+00, 6.33430e+00, 6.43030e+00, 6.63610e+00, 6.86260e+00, 7.11370e+00, 7.39330e+00, 7.54540e+00, 7.70680e+00, 8.06120e+00, 8.46520e+00, 8.93020e+00, 9.47130e+00, 1.01090e+01], [5.26880e+00, 5.31660e+00, 5.36580e+00, 5.46870e+00, 5.57820e+00, 5.69500e+00, 5.81980e+00, 5.88550e+00, 5.95540e+00, 6.09760e+00, 6.25300e+00, 6.42140e+00, 6.60440e+00, 6.70210e+00, 6.80430e+00, 7.02370e+00, 7.26550e+00, 7.53380e+00, 7.83310e+00, 7.99670e+00, 8.16940e+00, 8.55020e+00, 8.98510e+00, 9.48660e+00, 1.00713e+01, 1.07619e+01], [5.39660e+00, 5.44540e+00, 5.49570e+00, 5.60090e+00, 5.71280e+00, 5.83230e+00, 5.96010e+00, 6.02740e+00, 6.09720e+00, 6.24470e+00, 6.40410e+00, 6.57680e+00, 6.76470e+00, 6.86500e+00, 6.97000e+00, 7.19540e+00, 7.44420e+00, 7.72030e+00, 8.02860e+00, 8.19670e+00, 8.37530e+00, 8.76810e+00, 9.21810e+00, 9.73520e+00, 1.03399e+01, 1.10546e+01], [5.62550e+00, 5.67620e+00, 5.72840e+00, 5.83770e+00, 5.95410e+00, 6.07850e+00, 6.21160e+00, 6.28180e+00, 6.35460e+00, 6.50870e+00, 6.67520e+00, 6.85580e+00, 7.05260e+00, 7.15780e+00, 7.26790e+00, 7.50450e+00, 7.76600e+00, 8.05650e+00, 8.38130e+00, 8.55850e+00, 8.74690e+00, 9.16180e+00, 9.63670e+00, 1.01856e+01, 1.08270e+01, 1.15863e+01], [5.91700e+00, 5.97010e+00, 6.02480e+00, 6.13950e+00, 6.26180e+00, 6.39250e+00, 6.53270e+00, 6.60660e+00, 6.68330e+00, 6.84580e+00, 7.02180e+00, 7.21290e+00, 7.42130e+00, 7.53280e+00, 7.64960e+00, 7.90100e+00, 8.17910e+00, 8.48860e+00, 8.83520e+00, 9.02450e+00, 9.22600e+00, 9.67010e+00, 1.01793e+01, 1.07688e+01, 1.14590e+01, 1.22775e+01], [6.62100e+00, 6.68010e+00, 6.74110e+00, 6.86920e+00, 7.00620e+00, 7.15290e+00, 7.31070e+00, 7.39410e+00, 7.48070e+00, 7.66470e+00, 7.86440e+00, 8.08190e+00, 8.31890e+00, 8.44750e+00, 8.58140e+00, 8.87020e+00, 9.19080e+00, 9.54880e+00, 9.95110e+00, 1.01714e+01, 1.04062e+01, 1.09254e+01, 1.15229e+01, 1.22172e+01, 1.30332e+01, 1.40048e+01], [7.46200e+00, 7.52880e+00, 7.59770e+00, 7.74280e+00, 7.89820e+00, 8.06530e+00, 8.24540e+00, 8.34090e+00, 8.44020e+00, 8.65150e+00, 8.88160e+00, 9.13300e+00, 9.40900e+00, 9.55720e+00, 9.71320e+00, 1.00504e+01, 1.04259e+01, 1.08466e+01, 1.13211e+01, 1.15818e+01, 1.18601e+01, 1.24771e+01, 1.31898e+01, 1.40213e+01, 1.50024e+01, 1.61752e+01], [7.73620e+00, 7.80790e+00, 7.88210e+00, 8.03830e+00, 8.20610e+00, 8.38660e+00, 8.58160e+00, 8.68500e+00, 8.79270e+00, 9.02210e+00, 9.27240e+00, 9.54640e+00, 9.84770e+00, 1.00099e+01, 1.01805e+01, 1.05499e+01, 1.09622e+01, 1.14250e+01, 1.19480e+01, 1.22357e+01, 1.25432e+01, 1.32260e+01, 1.40164e+01, 1.49404e+01, 1.60330e+01, 1.73420e+01]]
# -> Problems for theta < 0.24 corresponding to int(Z) <= 10 ( theta_L_array(Z=11) = 0.24553 and theta_L_array(Z=10) = 0.23230 )
#    NEEDS SMOOTHING IN THIS CASE???
	B = maths_library.evaluate_value_from_2D_table(theta_LII_array,eta_LII_array,B_LII_table,theta,eta) 
	return max(0.,B);

def B_L(theta,eta):
# input  : theta = 1 - (I_R - I_0 )/I_NR where I_NR = theoretical non-relativistic ionization potential of L-shell electron(s)
#                                              I_R  = theoretical relativistic ionization potential of L-shell electron(s)
#                                              I_0  = experimental energy difference from the L-shell to the lowest-lying unoccupied state
#          eta   = me v^2 / (2 Zeff^2 Ry) where v is the projectile velocity and Zeff the effective charge seen by a LII-shell electron
# output : B     = L-shell contribution to the stopping number
#
	if (theta <= theta_L(2)):                                # i.e. if (int(Z) <= 2)
		n_LI = 0.
	elif (theta > theta_L(2)) and (theta <= theta_L(3)):     # i.e. if (int(Z) > 2) and (int(Z) <= 3)
		n_LI = 0.5
	else :                                                   # i.e. if (int(Z) > 3)
		n_LI = 1.
	if (theta <= theta_L(4)):                                # i.e. if (int(Z) <= 4)
		n_LII = 0.
	elif (theta > theta_L(4)) and (theta <= theta_L(9)):     # i.e. if (int(Z) > 4) and (int(Z) <= 9)
#
		if (theta > theta_L(4)) and (theta <= theta_L(5)):   # i.e. if int(Z) == 5
			n_LII = 0.5
		elif (theta > theta_L(5)) and (theta <= theta_L(6)): # i.e. if int(Z) == 6
			n_LII = 1.
		elif (theta > theta_L(6)) and (theta <= theta_L(7)): # i.e. if int(Z) == 7
			n_LII = 1.5
		elif (theta > theta_L(7)) and (theta <= theta_L(8)): # i.e. if int(Z) == 8
			n_LII = 2.
		elif (theta > theta_L(8)) and (theta <= theta_L(9)): # i.e. if int(Z) == 9
			n_LII = 2.5
	else :                                                   # i.e. if (int(Z) > 9)
		n_LII = 3.
	B = n_LI * B_LI(theta,eta) + ( n_LII * B_LII(theta,eta) )
	return B;

def eta_L(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : eta = (me/m) E / (Zeff^2 Ry)
#
#   Slater's rule :
	if (int(Z) <= 3):
		Zeff = float(Z)
	elif (int(Z) >4) and (int(Z) <=9):
		Zeff = float(Z) - ( 2. * 0.85 ) - ( float(Z-3) * 0.35 )
	else :
		Zeff = float(Z) - 4.15
	Zeff2 = Zeff * Zeff
	E     = nrj * 1.e6 / Ry
	a     = me / m
	eta   = a * E / Zeff2
	return eta;

def C_L(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C     = L-shell correction contribution to the stopping number
#
	theta    = theta_L(Z)
	S        = S_L(Z)
	T        = T_L(Z)
	U        = U_L(Z)
	eta      = eta_L(Z,m,nrj)
	ieta     = 1. / eta
	ieta_min = 0.1
	ieta_max = 0.6
	if (ieta < ieta_min) :
#   according to M. C. Walske, Phys. Rev. Vol. 88 No. 6 P. 1283 (1952):
#   C_L(Z,eta_K) = (U_L / eta_L) + (V_L / eta_L^2)
		ieta2 = ieta * ieta
		C = (U*ieta)
	elif (ieta >= ieta_min) and (ieta < ieta_max):
#   linear interpolation between the two expressions
		ieta_min2 = ieta_min * ieta_min
		Cmin      = (U*ieta_min)
		Bmax      = B_L(theta,1./ieta_max)
		Cmax      = ( S * math.log(1./ieta_max) ) + T - Bmax
		a         = (Cmax - Cmin) / (ieta_max - ieta_min)
		b         = Cmax - ( a * ieta_max )
		C         = (a * ieta) + b
	else :
#   according to M. S. Livingston and H. A. Bethe, Revs. Modern Phys. 9, 263 (1937) :
#   C_L(Z,eta_K) = S(Z) log(eta_L) + T(Z) - B_L(Z,eta_K)
		B = B_L(theta,eta)
		C = ( S * math.log(eta) ) + T - B
	return C;

###########
# model 1 #
###########

def V_M_model1(Z):
# input  : Z   = Material atomic number in ()
# output : parameter V_M - model 1 of Berger et al. 
#          Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
	if   (int(Z) <= 10 ):
		V = 0.
	elif (int(Z) > 10) and (int(Z) <= 32):
		V = (float(Z)-10.)/8.
	else:
		V = 18./8.
	return V;

def C_M_model1(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = M-shell correction contribution to the stopping number
#
	V_M      = V_M_model1(Z)
	z        = int(Z)-1
	s        = physical_constants.symbol
	H_M      = eval("physical_constants.M_shell_scale_parameter_model1_"+s[z])
	nrj_M    = H_M * nrj
	C_scaled = C_L(Z,m,nrj_M)
	C        = V_M * C_scaled 
	return C;

def V_N_model1(Z):
# input  : Z   = Material atomic number in ()
# output : parameter V_N - model 1 of Berger et al. 
#          Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
	if   (int(Z) <= 32 ):
		V = 0.
	elif (int(Z) > 32) and (int(Z) <= 60):
		V = (float(Z)-28.)/8.
	else:
		V = 32./8.
	return V;

def C_N_model1(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C     = N-shell correction contribution to the stopping number
#
	V_N       = V_N_model1(Z)
	z         = int(Z)-1
	s         = physical_constants.symbol
	H_N       = eval("physical_constants.N_shell_scale_parameter_model1_"+s[z])
	nrj_N     = H_N * nrj
	C_scaled  = C_L(Z,m,nrj_N)
	C         = V_N * C_scaled
	return C;

def V_OP_model1(Z):
# input  : Z   = Material atomic number in ()
# output : parameter V_OP - model 1 of Berger et al. 
#          Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
	if   (int(Z) <= 60 ):
		V = 0.
	else:
		V = (float(Z)-60.)/8.
	return V;

def H_OP_model1(Z):
# input  : Z   = Material atomic number in ()
# output : parameter H_OP - model 1 of Berger et al. 
#          Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
	if   (int(Z) <= 60 ):
		H = 0.
	else:
		H = 150.
	return H;

def C_OP_model1(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C     = OP-shells correction contribution to the stopping number
#
	V_OP     = V_OP_model1(Z)
	H_OP     = H_OP_model1(Z)
	nrj_OP   = H_OP * nrj
	C_scaled = C_L(Z,m,nrj)
	C        = V_OP * C_scaled
	return C;

###########
# model 2 #
###########

def theta_M_3s(Z):
# input  : Z     = Atomic number in ()
# output : theta = 9 I_0 / (Zeff^2 Ry) where I_0  = experimental K-shell electron ionization energy
#                                            Zeff = effective charge seen by K-shell electrons 
# according to H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983)
	if (int(Z) <= 12 ):
		theta = 0. 
	else :
		s        = physical_constants.symbol
		z        = int(Z) - 1
		Ioni     = physical_constants.M_shell_ionization_potential_model2 
		I_exp    = Ioni[0][z]
		# For Z = 26, 43, 54 and 61, A. Bearden and A. F. Burr, Rev. Modern Phys. Vol. 39 No.1 (1967)
		#                            don't provide the ionization energy.
		# Therefore, we extrapolate :  
		if (int(Z) >=17) and (I_exp == 0.) :
			x1 = Z - 1;y1 = Ioni[0][z-1]
			x2 = Z + 1;y2 = Ioni[0][z+1]
			I_exp =maths_library.linear_extrapolation_1D(x1,y1,x2,y2,Z)
		# Same thing for 10 < int(Z) < 17 :
		if (int(Z) > 10) and (I_exp == 0.) :
			x1 = 17;y1 = Ioni[0][16]
			x2 = 18;y2 = Ioni[0][17]
			I_exp = maths_library.linear_extrapolation_1D(x1,y1,x2,y2,Z)
			I_exp = max(0.,I_exp)
		Zeff_tab = eval("physical_constants.effective_charge_"+s[z])
		Zeff     = Zeff_tab[3]
		I_NR     = (Zeff**2.) * Ry / 9.
		a        = 1. 
		b        = 3. / 12.
		c        = ( (Zeff * alpha)**2. ) / 3.
		corr     = 1. + (  c * ( a - b ) )
		I_R      = I_NR * corr
		theta    = I_exp / I_R
	return theta;

def theta_M_3p(Z):
# input  : Z     = Atomic number in ()
# output : theta = 9 I_0 / (Zeff^2 Ry) where I_0  = experimental K-shell electron ionization energy
#                                            Zeff = effective charge seen by K-shell electrons
# according to H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983)
	if (int(Z) <= 18 ):
		theta = 0. 
	else :
		s          = physical_constants.symbol
		z          = int(Z) - 1
		Ioni       = physical_constants.M_shell_ionization_potential_model2 
		Zeff_tab   = eval("physical_constants.effective_charge_"+s[z])
		Zeff       = Zeff_tab[4]
		I_NR       = (Zeff**2.) * Ry / 9.
		# MII :
		I_exp_MII  = Ioni[1][z]
		# For 12 < int(Z) < 17, A. Bearden and A. F. Burr, Rev. Modern Phys. Vol. 39 No.1 (1967)
		#                  don't provide the ionization energy.
		# Therefore, we extrapolate :  
		if (int(Z) > 12) and (I_exp_MII == 0.) :
			x1 = 17;y1 = Ioni[1][16]
			x2 = 18;y2 = Ioni[1][17]
			I_exp_MII = maths_library.linear_extrapolation_1D(x1,y1,x2,y2,Z)
		a_MII      = 1. 
		b          = 3. / 12.
		c          = ( (Zeff * alpha)**2. ) / 3.
		corr_MII   = 1. + (  c * ( a_MII - b ) )
		I_MII      = I_NR * corr_MII
		n_MII      = 1. / 3.
		# MIII :
		I_exp_MIII = Ioni[2][z]
		# For 12 < int(Z) < 17, A. Bearden and A. F. Burr, Rev. Modern Phys. Vol. 39 No.1 (1967)
		#                  don't provide the ionization energy.
		# Therefore, we extrapolate :  
		if (int(Z) >=13) and (I_exp_MIII == 0.) :
			x1 = 17;y1 = Ioni[2][16]
			x2 = 18;y2 = Ioni[2][17]
			I_exp_MIII = maths_library.linear_extrapolation_1D(x1,y1,x2,y2,Z)
		a_MIII     = 0.5 
		corr_MIII  = 1. + (  c * ( a_MIII - b ) )
		I_MIII     = I_NR * corr_MIII
		n_MIII     = 2. / 3.
		# 3p :
		theta      = (n_MII * I_exp_MII / I_MII) + (n_MIII * I_exp_MIII / I_MIII)
	return theta;

def theta_M_3d(Z):
# input  : Z     = Atomic number in ()
# output : theta = 9 I_0 / (Zeff^2 Ry) where I_0  = experimental K-shell electron ionization energy
#                                            Zeff = effective charge seen by K-shell electrons 
# according to H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983)
	if (int(Z) <= 30 ):
		theta = 0. 
	else :
		s         = physical_constants.symbol
		z         = int(Z) - 1
		Ioni      = physical_constants.M_shell_ionization_potential_model2 
		Zeff_tab  = eval("physical_constants.effective_charge_"+s[z])
		Zeff      = Zeff_tab[6]
		I_NR      = (Zeff**2.) * Ry / 9.
		# MIV :
		I_exp_MIV = Ioni[3][z]
		# For Z = 54, A. Bearden and A. F. Burr, Rev. Modern Phys. Vol. 39 No.1 (1967)
		#             don't provide the ionization energy.
		# Therefore, we extrapolate :  
		if (int(Z) >=21) and (I_exp_MIV == 0.) :
			x1 = Z - 1;y1 = Ioni[3][z-1]
			x2 = Z + 1;y2 = Ioni[3][z+1]
			I_exp_MIV = maths_library.linear_extrapolation_1D(x1,y1,x2,y2,Z)
		a_MIV     = 0.5 
		b         = 3. / 12.
		c         = ( (Zeff * alpha)**2. ) / 3.
		corr_MIV  = 1. + (  c * ( a_MIV - b ) )
		I_MIV     = I_NR * corr_MIV
		n_MIV     = 2. / 5.
		# MV :
		I_exp_MV  = Ioni[4][z]
		a_MV      = 1. / 3. 
		corr_MV   = 1. + (  c * ( a_MV - b ) )
		I_MV      = I_NR * corr_MV
		n_MV      = 3. / 5.
		# 3d :
		theta     = (n_MIV * I_exp_MIV / I_MIV) + (n_MV * I_exp_MV / I_MV)
	return theta;

def S_M_3s(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - C_M(Z,eta_M)
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array   = [0.8100, 0.7200, 0.6300, 0.5400, 0.4500, 0.3600, 0.3150, 0.2700]
	S_array       = [1.3720, 1.4256, 1.4941, 1.5845, 1.7087, 1.8891, 2.0136, 2.1737]
	theta         = theta_M_3s(Z)
	S             = maths_library.evaluate_value_from_1D_table(theta,theta_array,S_array)
	return S; 

def T_M_3s(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - C_M(Z,eta_M)
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array   = [0.8100, 0.7200, 0.6300, 0.5400, 0.4500, 0.3600, 0.3150, 0.2700]
	T_array       = [5.0310, 5.2982, 5.6239, 6.0342, 6.5745, 7.3328, 7.8469, 8.5062]
	theta         = theta_M_3s(Z)
	T             = maths_library.evaluate_value_from_1D_table(theta,theta_array,T_array)
	return T; 

def U_M_3s(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - U(Z)/eta_M
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array = [0.8100, 0.7200, 0.6300, 0.5400, 0.4500, 0.3600, 0.3150, 0.2700]
	U_array     = [0.2548, 0.2552, 0.2553, 0.2551, 0.2547, 0.2541, 0.2539, 0.2539]
	theta       = theta_M_3s(Z)
	U           = maths_library.evaluate_value_from_1D_table(theta,theta_array,U_array)
	return U; 

def S_M_3p(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - C_M(Z,eta_M)
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array   = [0.81000, 0.72000, 0.63000, 0.54000, 0.45000, 0.36000, 0.31500]
	S_array       = [3.90580,  4.0897, 4.33490, 4.67270, 5.15870, 5.89950, 6.42800]
	theta         = theta_M_3p(Z)
	S             = maths_library.evaluate_value_from_1D_table(theta,theta_array,S_array)
	return S; 

def T_M_3p(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - C_M(Z,eta_M)
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array   = [0.81000, 0.72000, 0.63000, 0.54000,  0.45000,  0.36000,  0.31500]
	T_array       = [15.9256, 16.8309, 17.9452, 19.3686, 21.28250, 24.05640, 26.00200]
	theta         = theta_M_3p(Z)
	T             = maths_library.evaluate_value_from_1D_table(theta,theta_array,T_array)
	return T; 

def U_M_3p(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - U(Z)/eta_M
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array = [0.81000, 0.72000, 0.63000, 0.54000, 0.45000, 0.36000, 0.31500]
	U_array     = [0.64840,  0.6531, 0.65460, 0.65200, 0.64570, 0.64310, 0.65240]
	theta       = theta_M_3p(Z)
	U           = maths_library.evaluate_value_from_1D_table(theta,theta_array,U_array)
	return U; 

def S_M_3d(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - C_M(Z,eta_M)
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array   = [0.81000, 0.72000, 0.63000, 0.54000, 0.45000,  0.36000]
	S_array       = [5.82090, 6.09990, 6.52170, 7.19290, 8.33530, 10.46960]
	theta         = theta_M_3d(Z)
	S             = maths_library.evaluate_value_from_1D_table(theta,theta_array,S_array)
	return S; 

def T_M_3d(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - C_M(Z,eta_M)
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array   = [ 0.8100,  0.7200,  0.6300,  0.5400,  0.4500,  0.3600]
	T_array       = [28.0388, 29.8843, 32.2033, 35.2355, 39.4198, 45.6579]
	theta         = theta_M_3d(Z)
	T             = maths_library.evaluate_value_from_1D_table(theta,theta_array,T_array)
	return T; 

def U_M_3d(Z):
# input  : Z = Atomic number in ()
# output : S = S coefficient of asymptotic stopping number formula for the K-shell contribution to the stopping number
#              B_M(Z,eta_M) = S(Z) log(eta_M) + T(Z) - U(Z)/eta_M
#   according to Table I of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta_array = [0.81000, 0.72000, 0.63000, 0.54000, 0.45000, 0.36000]
	U_array     = [0.97910, 1.03100, 1.08480, 1.13120, 1.14320, 1.03970]
	theta       = theta_M_3d(Z)
	U           = maths_library.evaluate_value_from_1D_table(theta,theta_array,U_array)
	return U;

def C_M_3s_model2(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = M-shell (3s) electrons correction contribution to the stopping number
#                according to Tables I and II of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta           = theta_M_3s(Z)
	beta            = particle_beta_norm(m,nrj)
	beta2           = beta * beta
	s               = physical_constants.symbol
	z               = int(Z) - 1
	Zeff_tab        = eval("physical_constants.effective_charge_"+s[z])
	Zeff            = Zeff_tab[3]
	Zeff2           = Zeff * Zeff
	eta             = 18787. * beta2 / Zeff2
	if (eta < 2.5) :
		eta_array       =  [0.0100, 0.0126, 0.0158, 0.0200, 0.0251, 0.0316, 0.0398, 0.0501, 0.0631, 0.0794, 0.1000, 0.1259, 0.1585, 0.1995, 0.2512, 0.3162, 0.3981, 0.5012, 0.6310, 0.7943, 1.0000, 1.2589, 1.5849, 1.9953, 2.5119]
		theta_array     =  [ 0.810,  0.720,  0.630,  0.540,  0.450,  0.360,  0.315,  0.270]
		Bichsel_tableII = [[-1.291, -1.272, -1.267, -1.284, -1.337, -1.448, -1.535, -1.651], [-0.980, -0.953, -0.940, -0.946, -0.984, -1.072, -1.142, -1.235], [-0.678, -0.647, -0.629, -0.629, -0.656, -0.724, -0.780, -0.853], [-0.394, -0.363, -0.343, -0.340, -0.359, -0.413, -0.456, -0.511], [-0.138, -0.110, -0.092, -0.087, -0.102, -0.143, -0.175, -0.215], [ 0.081, 0.104, 0.199, 0.123, 0.112, 0.082, 0.058, 0.030], [ 0.256, 0.275, 0.287, 0.290, 0.281, 0.258, 0.241, 0.223], [ 0.387, 0.401, 0.410, 0.412, 0.405, 0.388, 0.377, 0.366], [ 0.477, 0.487, 0.494, 0.495, 0.489, 0.476, 0.469, 0.464], [ 0.530, 0.538, 0.543, 0.543, 0.538, 0.529, 0.524, 0.523], [ 0.554, 0.560, 0.563, 0.563, 0.559, 0.552, 0.549, 0.550], [ 0.554, 0.559, 0.561, 0.561, 0.558, 0.552, 0.550, 0.552], [ 0.537, 0.540, 0.542, 0.542, 0.539, 0.535, 0.534, 0.536], [ 0.507, 0.509, 0.511, 0.510, 0.508, 0.505, 0.504, 0.507], [ 0.468, 0.470, 0.471, 0.471, 0.469, 0.467, 0.466, 0.468], [ 0.425, 0.427, 0.427, 0.427, 0.426, 0.424, 0.423, 0.425], [ 0.380, 0.382, 0.382, 0.382, 0.381, 0.379, 0.379, 0.380], [ 0.336, 0.337, 0.337, 0.337, 0.336, 0.335, 0.335, 0.336], [ 0.294, 0.294, 0.294, 0.294, 0.294, 0.293, 0.292, 0.293], [ 0.253, 0.254, 0.254, 0.254, 0.253, 0.253, 0.252, 0.253], [ 0.216, 0.216, 0.217, 0.216, 0.216, 0.215, 0.215, 0.216], [ 0.182, 0.182, 0.182, 0.182, 0.182, 0.181, 0.181, 0.182], [ 0.151, 0.152, 0.152, 0.152, 0.151, 0.151, 0.151, 0.151], [ 0.124, 0.125, 0.125, 0.125, 0.124, 0.124, 0.124, 0.124], [ 0.101, 0.101, 0.101, 0.101, 0.101, 0.101, 0.101, 0.101]]
		C               = maths_library.evaluate_value_from_2D_table(theta_array,eta_array,Bichsel_tableII,theta,eta)
	else :
		C = U_M_3s(Z) / eta
	return C;

def C_M_3p_model2(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = M-shell (3p) electrons correction contribution to the stopping number
#                according to Tables I and III of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta           = theta_M_3p(Z)
	beta            = particle_beta_norm(m,nrj)
	beta2           = beta * beta
	s               = physical_constants.symbol
	z               = int(Z) - 1
	Zeff_tab        = eval("physical_constants.effective_charge_"+s[z])
	Zeff            = Zeff_tab[4]
	Zeff2           = Zeff * Zeff
	eta             = 18787. * beta2 / Zeff2
	if (eta < 2.5) :
		eta_array        =  [0.0100, 0.0126, 0.0158, 0.0200, 0.0251, 0.0316, 0.0398, 0.0501, 0.0631, 0.0794, 0.1000, 0.1259, 0.1585, 0.1995, 0.2512, 0.3162, 0.3981, 0.5012, 0.6310, 0.7943, 1.0000, 1.2589, 1.5849, 1.9953, 2.5119]
		theta_array      =  [  0.810,  0.720,  0.630,  0.540,  0.450,  0.360,  0.315]
		Bichsel_tableIII = [[-2.077, -2.029, -2.061, -2.225, -2.605, -3.341, -3.908], [-1.198, -1.118, -1.110, -1.221, -1.524, -2.139, -2.617], [-0.342, -0.239, -0.201, -0.269, -0.507, -1.016, -1.412], [ 0.471, 0.587, 0.644, 0.606, 0.420, 0.005, -0.318], [ 1.214, 1.333, 1.397, 1.379, 1.233, 0.898, 0.640], [ 1.857, 1.968, 2.032, 2.025, 1.910, 1.643, 1.442], [ 2.372, 2.470, 2.527, 2.525, 2.434, 2.222, 2.070], [ 2.740, 2.823, 2.871, 2.870, 2.797, 2.632, 2.521], [ 2.959, 3.027, 3.066, 3.064, 3.004, 2.877, 2.799], [ 3.040, 3.095, 3.125, 3.122, 3.072, 2.975, 2.924], [ 3.006, 3.050, 3.073, 3.069, 3.028, 2.954, 2.924], [ 2.884, 2.919, 2.937, 2.932, 2.897, 2.842, 2.827], [ 2.701, 2.728, 2.742, 2.736, 2.707, 2.667, 2.662], [ 2.475, 2.497, 2.508, 2.502, 2.478, 2.448, 2.450], [ 2.223, 2.240, 2.248, 2.243, 2.223, 2.200, 2.207], [ 1.955, 1.969, 1.975, 1.970, 1.954, 1.937, 1.946], [ 1.683, 1.695, 1.699, 1.695, 1.681, 1.669, 1.679], [ 1.419, 1.428, 1.432, 1.428, 1.417, 1.408, 1.417], [ 1.173, 1.180, 1.183, 1.180, 1.171, 1.164, 1.173], [ 0.952, 0.958, 0.960, 0.957, 0.950, 0.945, 0.953], [ 0.762, 0.766, 0.768, 0.766, 0.760, 0.756, 0.763], [ 0.603, 0.606, 0.608, 0.606, 0.601, 0.598, 0.604], [ 0.473, 0.476, 0.477, 0.476, 0.472, 0.470, 0.474], [ 0.370, 0.372, 0.373, 0.372, 0.369, 0.367, 0.371], [ 0.289, 0.291, 0.291, 0.290, 0.288, 0.286, 0.290]]
		C                = maths_library.evaluate_value_from_2D_table(theta_array,eta_array,Bichsel_tableIII,theta,eta)
	else :
		C = U_M_3p(Z) / eta
	return C;

def C_M_3d_model2(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = M-shell (3d) electrons correction contribution to the stopping number
#                according to Tables I and IV of H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983) :
	theta           = theta_M_3d(Z)
	beta            = particle_beta_norm(m,nrj)
	beta2           = beta * beta
	s               = physical_constants.symbol
	z               = int(Z) - 1
	Zeff_tab        = eval("physical_constants.effective_charge_"+s[z])
	Zeff            = Zeff_tab[6]
	Zeff2           = Zeff * Zeff
	eta             = 18787. * beta2 / Zeff2
	if (eta < 1.) :
		eta_array       =  [0.0100, 0.0126, 0.0158, 0.0200, 0.0251, 0.0316, 0.0398, 0.0501, 0.0631, 0.0794, 0.1000, 0.1259, 0.1585, 0.1995, 0.2512, 0.3162, 0.3981, 0.5012, 0.6310, 0.7943, 1.0000]
		theta_array     =  [      0.810,      0.720,      0.630,      0.540,      0.450,       0.360]
		Bichsel_tableIV = [[ 1.1870e+00, 1.7270e+00, 2.0740e+00, 1.0000e+00, 8.2400e-01, -2.8890e+00], [ 2.4850e+00, 3.0770e+00, 3.5040e+00, 9.7100e+02, 2.6040e+00, -6.8300e-01], [ 3.7550e+00, 4.3910e+00, 4.8870e+00, 3.5300e+00, 4.3010e+00, 1.4140e+00], [ 4.9760e+00, 5.6440e+00, 6.1920e+00, 5.0250e+00, 5.8760e+00, 3.3610e+00], [ 6.1170e+00, 6.8000e+00, 7.3790e+00, 6.4230e+00, 7.2810e+00, 5.1050e+00], [ 7.1310e+00, 7.8060e+00, 8.3940e+00, 7.6780e+00, 8.4600e+00, 6.5900e+00], [ 7.9580e+00, 8.6030e+00, 9.1750e+00, 8.7360e+00, 9.3550e+00, 7.7600e+00], [ 8.5320e+00, 9.1290e+00, 9.6670e+00, 9.5360e+00, 9.9160e+00, 8.5630e+00], [ 8.7990e+00, 9.3350e+00, 9.8250e+00, 1.0026e+01, 1.0109e+01, 8.9690e+00], [ 8.7300e+00, 9.1990e+00, 9.6350e+00, 1.0168e+01, 9.9270e+00, 8.9710e+00], [ 8.3340e+00, 8.7360e+00, 9.1150e+00, 9.9500e+00, 9.3960e+00, 8.6000e+00], [ 7.6590e+00, 7.9970e+00, 8.3210e+00, 9.3970e+00, 8.5800e+00, 7.9210e+00], [ 6.7850e+00, 7.0660e+00, 7.3380e+00, 8.5680e+00, 7.5700e+00, 7.0270e+00], [ 5.8080e+00, 6.0390e+00, 6.2650e+00, 7.5500e+00, 6.4680e+00, 6.0220e+00], [ 4.8220e+00, 5.0110e+00, 5.1970e+00, 6.4450e+00, 5.3700e+00, 5.0070e+00], [ 3.9020e+00, 4.0550e+00, 4.2080e+00, 5.3470e+00, 4.3540e+00, 4.0580e+00], [ 3.0950e+00, 3.2180e+00, 3.3420e+00, 4.3320e+00, 3.4640e+00, 3.2240e+00], [ 2.4200e+00, 2.5190e+00, 2.6200e+00, 3.4440e+00, 2.7200e+00, 2.5260e+00], [ 1.8770e+00, 1.9560e+00, 2.0370e+00, 2.7030e+00, 2.1190e+00, 1.9630e+00], [ 1.4500e+00, 1.5130e+00, 1.5780e+00, 2.1040e+00, 1.6450e+00, 1.5200e+00], [ 1.1200e+00, 1.1710e+00, 1.2220e+00, 1.6330e+00, 1.2760e+00, 1.1760e+00]]
		C               = maths_library.evaluate_value_from_2D_table(theta_array,eta_array,Bichsel_tableIV,theta,eta)
	else :
		C = U_M_3d(Z) / eta
	return C;

def C_M_model2(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = M-shell correction contribution to the stopping number
#
#   according to H. Bichsel, Phys. Rev. A, Vol. 28, No. 2 (1983)
	C = 0.
	if (int(Z) > 12):
		C = C + C_M_3s_model2(Z,m,nrj)
	if (int(Z) > 18):
		C = C + C_M_3p_model2(Z,m,nrj)
	if (int(Z) > 30):
		C = C + C_M_3d_model2(Z,m,nrj)
	return C;

def C_N1_model2(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = NI to NII subshell correction contribution to the stopping number
#
    # according to H. Bichsel "Stopping Power of Fast Charged Particles in Heavy Elements" 
    #              National Inst, of Standards and Technology, Report NIST IR-4550 (1991) -
    #              National Technical Information Service, Springfield, Virginia
	if (int(Z) >= 37):
		z     = int(Z) - 1
		J_MV  = physical_constants.M_shell_ionization_potential_model2[4][z]
		if (J_MV == 0):
			x1 =float(Z)-1.;y1 = physical_constants.M_shell_ionization_potential_model2[4][z-1]
			x2 =float(Z)+1.;y2 = physical_constants.M_shell_ionization_potential_model2[4][z+1]
			J_MV = maths_library.linear_interpolation_1D(x1,y1,x2,y2,float(Z))
		n_NI  = 2.
		J_NI  = physical_constants.N_shell_ionization_potential_model2[0][z]
		if (J_NI == 0.) :
			x1   = float(Z)-1.;y1 = physical_constants.N_shell_ionization_potential_model2[0][z-1]
			x2   = float(Z)+1.;y2 = physical_constants.N_shell_ionization_potential_model2[0][z+1]
			J_NI = maths_library.linear_interpolation_1D(x1,y1,x2,y2,float(Z))
		n_NII = 2.
		J_NII = physical_constants.N_shell_ionization_potential_model2[1][z]
		if (J_NII == 0.) :
			x1    = float(Z)-1.;y1 = physical_constants.N_shell_ionization_potential_model2[1][z-1]
			x2    = float(Z)+1.;y2 = physical_constants.N_shell_ionization_potential_model2[1][z+1]
			J_NII = maths_library.linear_interpolation_1D(x1,y1,x2,y2,float(Z))
		n_NIII = 4.	
		J_NIII = physical_constants.N_shell_ionization_potential_model2[2][z]
		if (J_NIII == 0.) :
			x1     = float(Z)-1.;y1 = physical_constants.N_shell_ionization_potential_model2[2][z-1]
			x2     = float(Z)+1.;y2 = physical_constants.N_shell_ionization_potential_model2[2][z+1]
			J_NIII = maths_library.linear_interpolation_1D(x1,y1,x2,y2,float(Z))	
		H1 = J_MV * n_NI / J_NI
		H1 = H1 + ( J_MV * n_NII  / J_NII )
		H1 = H1 + ( J_MV * n_NIII / J_NIII )
		H1 = H1 / ( n_NI + n_NII + n_NIII )
		V1 = 1.25
		C  = V1 * C_M_3d_model2(Z,m,H1*nrj)
	else :
		C = 0.
	return C;

def C_N2_model2(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = NIV to NV subshell correction contribution to the stopping number
#
    # according to H. Bichsel "Stopping Power of Fast Charged Particles in Heavy Elements" 
    #              National Inst, of Standards and Technology, Report NIST IR-4550 (1991) -
    #              National Technical Information Service, Springfield, Virginia
	if (int(Z) >= 41):
		z     = int(Z) - 1
		J_MV  = physical_constants.M_shell_ionization_potential_model2[4][z]
		if (J_MV == 0):
			x1 =float(Z)-1.;y1 = physical_constants.M_shell_ionization_potential_model2[4][z-1]
			x2 =float(Z)+1.;y2 = physical_constants.M_shell_ionization_potential_model2[4][z+1]
			J_MV = maths_library.linear_interpolation_1D(x1,y1,x2,y2,float(Z))
		n_NIV = 4.
		J_NIV = physical_constants.N_shell_ionization_potential_model2[3][z]
		if (J_NIV == 0.) :
			x1    = float(Z)-1.;y1 = physical_constants.N_shell_ionization_potential_model2[3][z-1]
			x2    = float(Z)+1.;y2 = physical_constants.N_shell_ionization_potential_model2[3][z+1]
			J_NIV = maths_library.linear_interpolation_1D(x1,y1,x2,y2,float(Z))
		n_NV = 6.
		J_NV = physical_constants.N_shell_ionization_potential_model2[4][z]
		if (J_NV == 0.) :
			x1    = float(Z)-1.;y1 = physical_constants.N_shell_ionization_potential_model2[4][z-1]
			x2    = float(Z)+1.;y2 = physical_constants.N_shell_ionization_potential_model2[4][z+1]
			J_NV  = maths_library.linear_interpolation_1D(x1,y1,x2,y2,float(Z))
		H2 = J_MV * n_NIV / J_NIV
		H2 = H2 + (J_MV * n_NV  / J_NV )
		H2 = H2 / (n_NIV + n_NV)
		V2 = 1.4
		C  = V2 * C_M_3p_model2(Z,m,H2*nrj)
	else :
		C = 0.
	return C;

def C_N3_model2(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : C   = all other shells correction contribution to the stopping number
#
    # according to H. Bichsel "Stopping Power of Fast Charged Particles in Heavy Elements" 
    #              National Inst, of Standards and Technology, Report NIST IR-4550 (1991) -
    #              National Technical Information Service, Springfield, Virginia
	if (int(Z) >= 46):
		if (int(Z) >= 73):
			V3 = (float(Z) - 46.) / 25.
			H3 = 13. 
		elif (int(Z) >= 60) and (int(Z) < 73):
			V3 = 2.3
			H3 = 25.
		else :
			V3 = 3.85
			H3 = 50.
		C  = V3 * C_M_3d_model2(Z,m,H3*nrj)
	else :
		C = 0.
	return C;

#########
# final #
#########

def Stopping_number_shell_correction(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction : C / Z in ()
#          according to Berger et al. - 
#          - Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
#
	C = 0.
	if (int(Z) > 2):
		C = C + max(0.,C_K(Z,m,nrj))
	if (int(Z) > 10):
		C = C + max(0.,C_L(Z,m,nrj))
	if (int(Z) < 64):
		if (int(Z) > 30):
			C = C + max(0.,C_M_model1(Z,m,nrj))
		if (int(Z) > 48):
			C = C + max(0.,C_N_model1(Z,m,nrj))
		if (int(Z) > 80):
			C = C + max(0.,C_OP_model1(Z,m,nrj))
	else:
		if (int(Z) > 30):
			C = C + max(0.,C_M_model2(Z,m,nrj))
		if (int(Z) > 48):
			C = C + max(0.,C_N1_model2(Z,m,nrj))
			C = C + max(0.,C_N2_model2(Z,m,nrj))
		if (int(Z) > 80):
			C = C + max(0.,C_N3_model2(Z,m,nrj))
	C = C / float(Z)
	return C;

######################################################################################
#                                 Low energy correction                              #
#                               C. Varelas and J. Biersack                           #
#  Reflection of Energetic Particles from Atomic or Ionic Chains in Single Crystals  #
#                           Nucl. Instr. Meth. 79 , 213 (1970)                       #
######################################################################################

def Stopping_power_low_energy_proton(material,nrj):
# input  : material = material name (string of characters)
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
#
	rho = get_density(material)
	A   = get_standard_atomic_weight(material)
	mi  = A * mu
	ni  = rho / mi
	[A1,A2,A3,A4,A5] = get_low_energy_proton_empirical_formula_A_parameters(material)
	T_s = nrj * 1.e3 / 1.0073
	# if (T_s < 10.) and (T_s >= 1.):
	# 	eps = A1 * (T_s**0.5)
	# else:
#	eps_low  = A2 * (T_s**0.45)
	eps_low  = A1 * (T_s**0.5)
	eps_high = (A3 / T_s) * math.log( 1. + (A4/T_s) + (A5*T_s) )
	eps      = eps_low * eps_high / (eps_low + eps_high)
	S        = eps * ni / 1.e21
	return S;

def Stopping_power_low_energy_alpha(material,nrj):
# input  : material = material name (string of characters)
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
#
	rho = get_density(material)
	A   = get_standard_atomic_weight(material)
	mi  = A * mu
	ni  = rho / mi
	[a1,a2,a3,a4,a5] = get_low_energy_alpha_empirical_formula_a_parameters(material)
	eps_low  = a1 * ((nrj*1.e3)**a2)
	eps_high = (a3 / nrj) * math.log( 1. + (a4/nrj) + (a5*nrj) )
	eps      = eps_low * eps_high / (eps_low + eps_high)
	S        = eps * ni / 1.e21
	return S;

def Stopping_power_low_energy(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
#
	rho = get_density(material)
	A   = get_standard_atomic_weight(material)
	mi  = A * mu
	ni  = rho / mi
	[a1,a2,a3,a4,a5] = get_low_energy_alpha_empirical_formula_a_parameters(material)
	eps = a1 * ((nrj*1.e3)**a2) # doesn t work
	S   = eps * ni / 1.e21
	return S;

######################################################################################
#                               Nuclear stopping power                               #
######################################################################################

def Scaled_Stopping_power_nuclear_Moliere(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Scaled stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993)
# Table 4.1
	S_nucl_Moliere = [[  1.0e+08,  8.0e+07,  6.0e+07,  5.0e+07,  4.0e+07,  3.0e+07,  2.0e+07,  1.5e+07,  1.0e+07,  8.0e+06,  6.0e+06,  5.0e+06,  4.0e+06,  3.0e+06,  2.0e+06,  1.5e+06,  1.0e+06,  8.0e+05,  6.0e+05,  5.0e+05,  4.0e+05,  3.0e+05,  2.0e+05,  1.5e+05,  1.0e+05,  8.0e+04,  6.0e+04,  5.0e+04,  4.0e+04,  3.0e+04,  2.0e+04,  1.5e+04,  1.0e+04,  8.0e+03,  6.0e+03,  5.0e+03,  4.0e+03,  3.0e+03,  2.0e+03,  1.5e+03,  1.0e+03,  8.0e+02,  6.0e+02,  5.0e+02,  4.0e+02,  3.0e+02,  2.0e+02,  1.5e+02,  1.0e+02,  8.0e+01,  6.0e+01,  5.0e+01,  4.0e+01,  3.0e+01,  1.5e+01,  1.0e+01,  8.0e+00,  6.0e+00,  5.0e+00,  4.0e+00,  3.0e+00,  2.0e+00,  1.5e+00,  1.0e+00,  8.0e-01,  6.0e-01,  5.0e-01,  4.0e-01,  3.0e-01,  2.0e-01,  1.5e-01,  1.0e-01,  8.0e-02,  6.0e-02,  5.0e-02,  4.0e-02,  3.0e-02,  2.0e-02,  1.5e-02,  1.0e-02,  8.0e-03,  6.0e-03,  5.0e-03,  4.0e-03,  3.0e-03,  2.0e-03,  1.5e-03,  1.0e-03,  8.0e-04,  6.0e-04,  5.0e-04,  4.0e-04,  3.0e-04,  2.0e-04,  1.5e-04,  1.0e-04,  8.0e-05,  6.0e-05,  5.0e-05,  4.0e-05,  3.0e-05,  2.0e-05,  1.5e-05,  1.0e-05],
	                  [5.831e-08,7.288e-08,9.719e-08,1.166e-07,1.457e-07,1.942e-07,2.916e-07,3.887e-07,5.833e-07,7.287e-07,9.712e-07,1.166e-06,1.457e-06,1.941e-06,2.911e-06,3.878e-06,5.810e-06,7.262e-06,9.663e-06,1.157e-05,1.442e-05,1.913e-05,2.845e-05,3.762e-05,5.554e-05,6.866e-05,9.020e-05,1.070e-04,1.319e-04,1.722e-04,2.499e-04,3.248e-04,4.688e-04,5.729e-04,7.411e-04,8.718e-04,1.063e-03,1.370e-03,1.955e-03,2.511e-03,3.563e-03,4.314e-03,5.511e-03,6.430e-03,7.756e-03,9.855e-03,1.375e-02,1.736e-02,2.395e-02,2.850e-02,3.552e-02,4.073e-02,4.802e-02,5.904e-02,9.426e-02,1.210e-01,1.377e-01,1.611e-01,1.768e-01,1.968e-01,2.235e-01,2.613e-01,2.871e-01,3.199e-01,3.354e-01,3.523e-01,3.609e-01,3.693e-01,3.766e-01,3.803e-01,3.788e-01,3.711e-01,3.644e-01,3.530e-01,3.444e-01,3.323e-01,3.144e-01,2.854e-01,2.629e-01,2.298e-01,2.115e-01,1.883e-01,1.741e-01,1.574e-01,1.372e-01,1.116e-01,9.559e-02,7.601e-02,6.668e-02,5.605e-02,5.008e-02,4.352e-02,3.617e-02,2.768e-02,2.279e-02,1.723e-02,1.473e-02,1.200e-02,1.052e-02,8.950e-03,7.246e-03,5.358e-03,4.313e-03,3.166e-03]]
	A = get_standard_atomic_weight(material)
	Z = get_atomic_number(material) 
	M = A * mu
	a = float(Z)**(1./3.)
	T = 32.536e3 * nrj / ( a * (1. + (m/M)) * float(Z) )
	S = maths_library.evaluate_value_from_1D_table(T,S_nucl_Moliere[0],S_nucl_Moliere[1])
	return S;

def Stopping_power_nuclear_Moliere(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993)
# Eq. (4.15)
	rho    = get_density(material)
	A      = get_standard_atomic_weight(material)
	Z      = get_atomic_number(material)
	M      = A * mu
	a      = float(Z)**(1./3.)
	S_nucl = Scaled_Stopping_power_nuclear_Moliere(material,z,m,nrj)
	S      = rho * 5105.3 * float(Z) * S_nucl / ( a * ( 1. + (M/m) ) * A )
	return S;

def Scaled_Stopping_power_nuclear_Ziegler(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Scaled stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993)
# Eq. (4.14) and (4.16)
	A = get_standard_atomic_weight(material)
	Z = get_atomic_number(material)
	M = A * mu
	a = (float(z)**0.23) +(float(Z)**0.23)
	T = 32.536e3 * nrj / ( a * (1. + (m/M)) * float(Z) * float(z) )
	if (T < 30):
		Num = 0.5 * math.log(1.+(1.1383*T))
		Den = T + (0.01321*(T**0.21226)) + (0.1959*(T**0.5))
		S   = Num / Den
	else :
		S = 0.5 * math.log(T) / T
	return S;

def Stopping_power_nuclear_Ziegler(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993)
# Eq. (4.15)
	rho    = get_density(material)
	A      = get_standard_atomic_weight(material)
	Z      = get_atomic_number(material)
	M      = A * mu
	a      = (float(z)**0.23) +(float(Z)**0.23)
	S_nucl = Scaled_Stopping_power_nuclear_Ziegler(material,z,m,nrj)
	S      = rho * S_nucl * 5105.3 * float(z) * float(Z) / ( a * ( 1. + (M/m) ) * A )
	return S;

######################################################################################
#                      Compounds - Bragg's additivity rule                           #
#     Berger et al. - Stopping powers and ranges for protons and alpha particles -   #
#                                 ICRU report 49 (1993)                              #
######################################################################################

# Barkas correction for elements in compounds :

def scaled_minimum_impact_parameter_for_elements_in_compounds(Z):
# input  : Z = Compound element atomic number in ()
# output : b = Scaled minimum impact parameter to be used in
#              the Barkas correction on the stopping number
#   according to the Table 2.3 of Berger et al. 
#   - Stopping powers and ranges for protons and alpha particles -
#   ICRU report 49 (1993)
	if (int(Z) == 1):
		b = 1.8
	elif (int(Z) == 2):
		b = 0.6
	elif (int(Z) > 2) and (int(Z) <= 10) :
		b = 1.8
	elif (int(Z) > 10) and (int(Z) <= 17):
		b = 1.4
	elif (int(Z) == 18):
		b = 1.8
	elif (int(Z) > 18) and (int(Z) <= 25):
		b = 1.4
	elif (int(Z) > 25) and (int(Z) <= 50):
		b = 1.35
	else :
		b = 1.3
	return b;

def Barkas_correction_for_elements_in_compounds(Z,m,nrj):
# input  : Z   = Material atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction : 
#          L1(beta) = g F_ARB( b/sqrt(x) ) / ( sqrt(Z) * x^(3/2) ) in ()
#          where g = constant = 1.29 and x = (beta/alpha)^2 / Z with alpha the 
#          fine structure constant according to Eq (2.8) of Berger et al. 
#          - Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993)
	g     = 1.29
	beta  = particle_beta_norm(m,nrj)
	s     = physical_constants.symbol
	b     = scaled_minimum_impact_parameter_for_elements_in_compounds(Z)
	x     = ( ( beta / alpha )**2. ) / float(Z)
	w = b / (x**0.5)
	F_ARB = maths_library.Ashley_function(w)
	L1    = g * F_ARB / ( ( Z * (x**3.) )**0.5 )
	return L1;

def Stopping_number_Barkas_correction_for_elements_in_compounds(Z,z,m,nrj):
# input  : Z   = Material atomic number in ()
#          z   = projectile atomic number in ()
#          m   = projectile mass in g
#          nrj = projectile energy in (MeV)
# output : Stopping number correction : L1(beta) in ()
#          according to Berger et al. - 
#          - Stopping powers and ranges for protons and alpha particles -
#          ICRU report 49 (1993) 
	if (int(Z) < 64) and int(Z) != 47:
		L1 = Barkas_correction_for_elements_in_compounds(Z,m,nrj)
	elif int(Z) == 47:
		L1 = Bichsel_correction_silver(m,nrj)
	else :
		L1 = Bichsel_correction(Z,m,nrj)
	L1 = float(z) * L1 
	return L1;

# Bethe stopping power for elements in compounds :

def mean_excitation_energy_for_elements_in_compounds(Z):
# input  : Z = Compound element atomic number in ()
# output : I = Mean excitation energy in (eV) to be used in 
#          the Bethe stopping number
#   according to the Table 2.11 of Berger et al. 
#   - Stopping powers and ranges for protons and alpha particles -
#   ICRU report 49 (1993)
	if (int(Z) == 1):
		I = 19.2
	elif (int(Z) == 6):
		I = 81.
	elif (int(Z) == 7):
		I = 82.
	elif (int(Z) == 8):
		I = 106.
	elif (int(Z) == 9):
		I = 112.
	elif (int(Z) == 17):
		I = 180.
	else :
		mat = eval("physical_constants.name_"+s[int(Z)-1])
		I = get_mean_excitation_energy(mat)
		I = 1.13 * I
	return I;

def Stopping_number_Bethe_for_elements_in_compounds(material,m,nrj):
# input  : material = material name (string of characters)
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : Bethe stopping number L0(beta) in () :
	mec2  = me * c * c
	beta  = particle_beta_norm(m,nrj)
	beta2 = beta * beta
	x     = 1. - beta2
	Wm    = largest_possible_energy_loss_in_single_collision(m,nrj)
	Z     = get_atomic_number(material)
	I     = mean_excitation_energy_for_elements_in_compounds(Z) * eV
	L0    = ( 0.5 * math.log( 2. * mec2 * Wm *beta2 / x ) ) - math.log(I) - beta2
	return L0;

def Stopping_power_Bethe_for_elements_in_compounds(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
	# Bethe stopping number
	L0 = Stopping_number_Bethe_for_elements_in_compounds(material,m,nrj)
	Z  = get_atomic_number(material)
	# Shell correction
	C  = Stopping_number_shell_correction(Z,m,nrj)
	# We neglect density effects assuming non relativistic projectile energy
	D = 0.
	L0 = L0 - C - D
	# Barkas correction
	L1 = Stopping_number_Barkas_correction_for_elements_in_compounds(Z,z,m,nrj)
	# Bloch correction
	L2 = Stopping_number_Bloch_correction(z,m,nrj)
	#
	L  = L0 + L1 + L2
	# electron density
	rho = get_density(material)
	A   = get_standard_atomic_weight(material)
	mi  = A * mu
	ni  = rho / mi
	ne  = float(Z) * ni
	# factor
	v     = particle_velocity_norm(m, nrj)
	v2    = v * v
	qe2   = qe * qe
	qe4   = qe2 * qe2
	z2    = float(z) * float(z)
	alpha = 4. * pi * ne * qe4 * z2 / (me * v2)
	#
	S = alpha * L / MeV
	return S;

# Stopping power in compounds :

def Stopping_power_for_elements_in_compounds(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
	if (int(z) == 1) :
		[T1,T2] = get_low_energy_proton_empirical_formula_T_parameters(material)
		Z_mat = get_atomic_number(material)
		if (Z_mat != 29) and (Z_mat != 32) and (Z_mat != 36) and (Z_mat != 47) and (Z_mat != 51) and (Z_mat != 79):
			T1 = T2 
		if (nrj <= T1):
			S = Stopping_power_low_energy_proton(material, nrj)
		elif (nrj > T2):
			S = Stopping_power_Bethe_for_elements_in_compounds(material,z,m,nrj)
		else :
			S1    = Stopping_power_low_energy_proton(material, T1)
			S2    = Stopping_power_Bethe_for_elements_in_compounds(material,z,m,T2)
			S     = maths_library.linear_interpolation_1D(T1,S1,T2,S2,nrj)
		S = S + Stopping_power_nuclear_Moliere(material,z,m,nrj)
	elif (int(z) == 2) :
		[T1,T2] = get_low_energy_alpha_empirical_formula_T_parameters(material)
		if (nrj <= T1):
			S = Stopping_power_low_energy_alpha(material, nrj)
		elif (nrj > T2):
			S = Stopping_power_Bethe_for_elements_in_compounds(material,z,m,nrj)
		else :
			beta  = particle_beta_norm(m,nrj)
			beta2 = beta * beta
			S1    = beta2 * Stopping_power_low_energy_alpha(material, T1)
			S2    = beta2 * Stopping_power_Bethe_for_elements_in_compounds(material,z,m,T2)
			S     = maths_library.linear_interpolation_1D(T1,S1,T2,S2,nrj)/beta2
		S = S + Stopping_power_nuclear_Ziegler(material,z,m,nrj)
	else :
		# if (nrj <= 1.e1) :
		# 	S = Stopping_power_low_energy(Z,z,m,nrj)
		# else :
		# 	S = Stopping_power_Bethe(Z,z,m,nrj)
		S = Stopping_power_Bethe_for_elements_in_compounds(material,z,m,nrj)
		S = S + Stopping_power_nuclear_Ziegler(material,z,m,nrj) 
	return S;

def Stopping_power_in_compounds(material,z,m,nrj):
# input  : material = material name (string of characters)
#          z        = projectile atomic number in ()
#          m        = projectile mass in g
#          nrj      = projectile energy in (MeV)
# output : S        = Stopping power in (MeV/cm) of a particle projectile with mass m , electrical charge z 
#                 and kinetic energy nrj << m c^2 propagating in a material of atomic number Z
	if (material == 'GafChromic HD-V2 (solid at room temperature)'):
		rho = 1.2 # density in g/cm3
		# Hydrogen : n_X % of atoms X with atomic number Z_X and atomic weight A_X
		n_H    = 58.4
		Z_H    = 1
		A_H    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_H)-1])
		rho_H  = eval("physical_constants.density_"+s[int(Z_H)-1])
		w_tot  = n_H * A_H
		mat_H  = eval("physical_constants.name_"+s[int(Z_H)-1])
		S_H    = Stopping_power_for_elements_in_compounds(mat_H,z,m,nrj)
		# Lithium
		n_Li   = 0.6
		Z_Li   = 3
		A_Li   = eval("physical_constants.standard_atomic_weight_"+s[int(Z_Li)-1])
		rho_Li = eval("physical_constants.density_"+s[int(Z_Li)-1])
		w_tot  = w_tot + (n_Li * A_Li)
		mat_Li = eval("physical_constants.name_"+s[int(Z_Li)-1])
		S_Li   = Stopping_power_for_elements_in_compounds(mat_Li,z,m,nrj)
		# Carbon
		n_C    = 27.9
		Z_C    = 6
		A_C    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_C)-1])
		rho_C  = eval("physical_constants.density_"+s[int(Z_C)-1])[0] # amorpheous Carbon
		w_tot  = w_tot + n_C * A_C
		mat_C  = eval("physical_constants.name_"+s[int(Z_C)-1])[0]    # amorpheous Carbon
		S_C    = Stopping_power_for_elements_in_compounds(mat_C,z,m,nrj)
		# Nitrogen
		n_N    = 0.1
		Z_N    = 7
		A_N    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_N)-1])
		rho_N  = eval("physical_constants.density_"+s[int(Z_N)-1])
		w_tot  = w_tot + (n_N * A_N)
		mat_N  = eval("physical_constants.name_"+s[int(Z_N)-1])
		S_N    = Stopping_power_for_elements_in_compounds(mat_N,z,m,nrj)
		# Oxygen
		n_O    = 11.7
		Z_O    = 8
		A_O    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_O)-1])
		rho_O  = eval("physical_constants.density_"+s[int(Z_O)-1])
		w_tot  = w_tot + (n_O * A_O)
		mat_O  = eval("physical_constants.name_"+s[int(Z_O)-1])
		S_O    = Stopping_power_for_elements_in_compounds(mat_O,z,m,nrj)
		# Sodium
		n_Na   = 0.5
		Z_Na   = 11
		A_Na   = eval("physical_constants.standard_atomic_weight_"+s[int(Z_Na)-1])
		rho_Na = eval("physical_constants.density_"+s[int(Z_Na)-1])
		w_tot  = w_tot + (n_Na * A_Na)
		mat_Na = eval("physical_constants.name_"+s[int(Z_Na)-1])
		S_Na   = Stopping_power_for_elements_in_compounds(mat_Na,z,m,nrj)
		# Aluminum
		n_Al   = 0.3
		Z_Al   = 13
		A_Al   = eval("physical_constants.standard_atomic_weight_"+s[int(Z_Al)-1])
		rho_Al = eval("physical_constants.density_"+s[int(Z_Al)-1])
		w_tot  = w_tot + (n_Al * A_Al)
		mat_Al = eval("physical_constants.name_"+s[int(Z_Al)-1])
		S_Al   = Stopping_power_for_elements_in_compounds(mat_Al,z,m,nrj)
		# Chlorine
		n_Cl   = 0.6
		Z_Cl   = 17
		A_Cl   = eval("physical_constants.standard_atomic_weight_"+s[int(Z_Cl)-1])
		rho_Cl = eval("physical_constants.density_"+s[int(Z_Cl)-1])
		w_tot  = w_tot + (n_Cl * A_Cl)
		mat_Cl = eval("physical_constants.name_"+s[int(Z_Cl)-1])
		S_Cl   = Stopping_power_for_elements_in_compounds(mat_Cl,z,m,nrj)
		# Weights :
		w_H  = n_H  * A_H  * rho / ( w_tot * rho_H  )
		w_Li = n_Li * A_Li * rho / ( w_tot * rho_Li )
		w_C  = n_C  * A_C  * rho / ( w_tot * rho_C  )
		w_N  = n_N  * A_N  * rho / ( w_tot * rho_N  )
		w_O  = n_O  * A_O  * rho / ( w_tot * rho_O  )
		w_Na = n_Na * A_Na * rho / ( w_tot * rho_Na )
		w_Al = n_Al * A_Al * rho / ( w_tot * rho_Al )
		w_Cl = n_Cl * A_Cl * rho / ( w_tot * rho_Cl )
		# Stopping power
		S = w_H  * S_H
		S = S + ( w_Li * S_Li )
		S = S + ( w_C  * S_C  )
		S = S + ( w_N  * S_N  )
		S = S + ( w_O  * S_O  )
		S = S + ( w_Na * S_Na )
		S = S + ( w_Al * S_Al )
		S = S + ( w_Cl * S_Cl )
	elif (material == 'Polyester (solid at room temperature)'):
		rho = 1.35 # density in g/cm3
		# Hydrogen : n_X % of elements X with atomic number Z_X and atomic weight A_X
		n_H    = 36.4
		Z_H    = 1
		A_H    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_H)-1])
		rho_H  = eval("physical_constants.density_"+s[int(Z_H)-1])
		w_tot  = n_H * A_H
		mat_H  = eval("physical_constants.name_"+s[int(Z_H)-1])
		S_H    = Stopping_power_for_elements_in_compounds(mat_H,z,m,nrj)
		# Carbon
		n_C    = 45.5
		Z_C    = 6
		A_C    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_C)-1])
		rho_C  = eval("physical_constants.density_"+s[int(Z_C)-1])[0] # amorpheous Carbon
		w_tot  = w_tot + n_C * A_C
		mat_C  = eval("physical_constants.name_"+s[int(Z_C)-1])[0]    # amorpheous Carbon
		S_C    = Stopping_power_for_elements_in_compounds(mat_C,z,m,nrj)
		# Oxygen
		n_O    = 18.2
		Z_O    = 8
		A_O    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_O)-1])
		rho_O  = eval("physical_constants.density_"+s[int(Z_O)-1])
		w_tot  = w_tot + (n_O * A_O)
		mat_O  = eval("physical_constants.name_"+s[int(Z_O)-1])
		S_O    = Stopping_power_for_elements_in_compounds(mat_O,z,m,nrj)
		# Weights :
		w_H  = n_H  * A_H  * rho / ( w_tot * rho_H  )
		w_C  = n_C  * A_C  * rho / ( w_tot * rho_C  )
		w_O  = n_O  * A_O  * rho / ( w_tot * rho_O  )
		# Stopping power
		S = w_H  * S_H
		S = S + ( w_C  * S_C  )
		S = S + ( w_O  * S_O  )
	elif (material == 'Plastic scintillator (vinyltoluene-based solid at room temperature)'):
		rho    = 1.032 # density in g/cm3
		#
		den_H  = 5.23e22 # H atoms / cm^3
		den_C  = 4.74e22 # C atoms / cm^3
		# Hydrogen : n_X % of elements X with atomic number Z_X and atomic weight A_X
		n_H    = den_H / ( den_C + den_H )
		Z_H    = 1
		A_H    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_H)-1])
		m_H    = A_H * mu
		rho_H  = 5.23e22 * m_H
		w_tot  = n_H * A_H
		mat_H  = eval("physical_constants.name_"+s[int(Z_H)-1])
		S_H    = Stopping_power_for_elements_in_compounds(mat_H,z,m,nrj)
		# Carbon
		n_C    = den_C / ( den_C + den_H )
		Z_C    = 6
		A_C    = eval("physical_constants.standard_atomic_weight_"+s[int(Z_C)-1])
		m_C    = A_C * mu
		rho_C  = eval("physical_constants.density_"+s[int(Z_C)-1])[0] # amorpheous Carbon
		w_tot  = n_C * A_C
		mat_C  = eval("physical_constants.name_"+s[int(Z_C)-1])[0]    # amorpheous Carbon
		S_C    = Stopping_power_for_elements_in_compounds(mat_C,z,m,nrj)
		# Weights :
		w_H  = n_H  * A_H  * rho / ( w_tot * rho_H  )
		w_C  = n_C  * A_C  * rho / ( w_tot * rho_C  )
		# Stopping power
		S = w_H  * S_H
		S = S + ( w_C  * S_C  )
	elif (material == 'Water (liquid at room temperature)'):
		rho    = 1. # density in g/cm3
		# Hydrogen :
		n_H    = 0.111894
		Z_H    = 1
		rho_H  = eval("physical_constants.density_"+s[int(Z_H)-1])
		mat_H  = eval("physical_constants.name_"+s[int(Z_H)-1])
		S_H    = Stopping_power_for_elements_in_compounds(mat_H,z,m,nrj)
		# Oxygen
		n_O    = 0.888106
		Z_O    = 8
		rho_O  = eval("physical_constants.density_"+s[int(Z_O)-1])
		mat_O  = eval("physical_constants.name_"+s[int(Z_O)-1])
		S_O    = Stopping_power_for_elements_in_compounds(mat_O,z,m,nrj)
		# Weights :
		w_H  = n_H  * rho / rho_H
		w_O  = n_O  * rho / rho_O
		# Stopping power
		S = w_H  * S_H
		S = S + ( w_O  * S_O )
	elif (material == 'Air (dry gas at room temperature)') :
		rho    = 1.20484e-3 # density in g/cm3
		# Carbon :
		n_C    = 0.000124
		Z_C    = 6
		rho_C  = eval("physical_constants.density_"+s[int(Z_C)-1])[0]
		mat_C  = eval("physical_constants.name_"+s[int(Z_C)-1])[0]
		S_C    = Stopping_power_for_elements_in_compounds(mat_C,z,m,nrj)
		# Nitrogen
		n_N    = 0.755267
		Z_N    = 7
		rho_N  = eval("physical_constants.density_"+s[int(Z_N)-1])
		mat_N  = eval("physical_constants.name_"+s[int(Z_N)-1])
		S_N    = Stopping_power_for_elements_in_compounds(mat_N,z,m,nrj)
		# Oxygen
		n_O    = 0.231781
		Z_O    = 8
		rho_O  = eval("physical_constants.density_"+s[int(Z_O)-1])
		mat_O  = eval("physical_constants.name_"+s[int(Z_O)-1])
		S_O    = Stopping_power_for_elements_in_compounds(mat_O,z,m,nrj)
		# Argon
		n_Ar   = 0.012827
		Z_Ar   = 18
		rho_Ar = eval("physical_constants.density_"+s[int(Z_Ar)-1])
		mat_Ar = eval("physical_constants.name_"+s[int(Z_Ar)-1])
		S_Ar   = Stopping_power_for_elements_in_compounds(mat_Ar,z,m,nrj)
		# Weights :
		w_C  = n_C  * rho / rho_C
		w_N  = n_N  * rho / rho_N
		w_O  = n_O  * rho / rho_O
		w_Ar = n_Ar * rho / rho_Ar
		# Stopping power
		S = w_C  * S_C
		S = S + ( w_N  * S_N  )
		S = S + ( w_O  * S_O  )
		S = S + ( w_Ar * S_Ar )
	elif (material == 'Glass SiO2 (solid at room temperature)') :
		rho    = 2.32 # density in g/cm3
		# Oxygen
		n_O    = 0.532565
		Z_O    = 8
		rho_O  = eval("physical_constants.density_"+s[int(Z_O)-1])
		mat_O  = eval("physical_constants.name_"+s[int(Z_O)-1])
		S_O    = Stopping_power_for_elements_in_compounds(mat_O,z,m,nrj)
		# Silicon : 
		n_Si   = 0.467435
		Z_Si   = 14
		rho_Si = eval("physical_constants.density_"+s[int(Z_Si)-1])
		mat_Si = eval("physical_constants.name_"+s[int(Z_Si)-1])
		S_Si   = Stopping_power_for_elements_in_compounds(mat_Si,z,m,nrj)
		# Weights :
		w_O  = n_O  * rho / rho_O
		w_Si = n_Si * rho / rho_Si
		# Stopping power
		S = w_O  * S_O
		S = S + ( w_Si  * S_Si )		
	else :
		S = 0
		print('Unknown compound material')
	return S;
