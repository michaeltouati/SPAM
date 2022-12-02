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

######################################################################################
#                          Physical_constants constants
######################################################################################

# Atomic mass unit [g] :
atomic_mass_unit          = 1.6605e-24      
# Proton mass [g] :                     
proton_mass               = 1.6726e-24           
# Electron mass [g] :                
electron_mass             = 9.1094e-28
# Elementary charge [statcoulomb] :                  
elementary_charge         = 4.8032e-10
# Speed of light in vacuum [cm/s] :                           
light_speed_in_vacuum     = 2.9979e10
# mp c^2 en MeV :                             
proton_mass_energy        = 938.272088
# Boltzmann constant [erg/K] :  
Boltzmann_constant        = 1.3807e-16
# Planck constant [erg.s] :                           
Planck_constant           = 6.6261e-27                           
# Fine Structure constant :
fine_structure_constant   = 7.2972e-3    
# classical electron radius [cm]
classical_electron_radius = 2.8179e-13
# Rydberg constant in [eV]
Rydberg_constant          = 13.605693123

######################################################################################
#                              Conversion factors
######################################################################################

# 1 J in (erg) :
Joules_in_erg = 1.e7
# 1 eV = in (erg) :
eV_in_erg     = 1.6022e-19 * Joules_in_erg
# 1 keV in (erg) :             
keV_in_erg    = 1.e3 * eV_in_erg
# 1 MeV in (erg) :             
MeV_in_erg    = 1.e6 * eV_in_erg

# 1 fs in s :                
fs = 1.e-15
# 1 micron in cm :
microns = 1.e-4
# Ambiant temperature [K] :
Tamb=300

######################################################################################
#                              Material properties
######################################################################################

#################
# Atomic symbol #
#################
# create an array symbol providing the corresponding atomic symbol : for example, 
# if the atomic number Z is Z=1, then symbol[Z-1] provides the atomic symbol "H" of Hydrogen. 
symbol=[]
symbol.append("H")
symbol.append("He") 
symbol.append("Li")
symbol.append("Be")
symbol.append("B")
symbol.append("C")
symbol.append("N")
symbol.append("O")
symbol.append("F")
symbol.append("Ne")
symbol.append("Na")
symbol.append("Mg")
symbol.append("Al")
symbol.append("Si")
symbol.append("P")
symbol.append("S")
symbol.append("Cl")
symbol.append("Ar")
symbol.append("K")
symbol.append("Ca")
symbol.append("Sc")
symbol.append("Ti")
symbol.append("V")
symbol.append("Cr")
symbol.append("Mn")
symbol.append("Fe")
symbol.append("Co")
symbol.append("Ni")
symbol.append("Cu")
symbol.append("Zn")
symbol.append("Ga")
symbol.append("Ge")
symbol.append("As")
symbol.append("Se")
symbol.append("Br")
symbol.append("Kr")
symbol.append("Rb")
symbol.append("Sr")
symbol.append("Y")
symbol.append("Zr")
symbol.append("Nb")
symbol.append("Mo")
symbol.append("Tc") 
symbol.append("Ru")
symbol.append("Rh")
symbol.append("Pd")
symbol.append("Ag")
symbol.append("Cd")
symbol.append("In")
symbol.append("Sn")
symbol.append("Sb")
symbol.append("Te")
symbol.append("I")
symbol.append("Xe")
symbol.append("Cs")
symbol.append("Ba")
symbol.append("La")
symbol.append("Ce")
symbol.append("Pr")
symbol.append("Nd")
symbol.append("Pm")
symbol.append("Sm")
symbol.append("Eu")
symbol.append("Gd")
symbol.append("Tb")
symbol.append("Dy")
symbol.append("Ho")
symbol.append("Er")
symbol.append("Tm")
symbol.append("Yb")
symbol.append("Lu")
symbol.append("Hf")
symbol.append("Ta")
symbol.append("W")
symbol.append("Re")
symbol.append("Os")
symbol.append("Ir")
symbol.append("Pt")
symbol.append("Au")
symbol.append("Hg")
symbol.append("Tl")
symbol.append("Pb")

################
# Element name #
################
name_H  = "Hydrogen (gas H2 at room temperature)"
name_He = "Helium (noble gas at room temperature)"
name_Li = "Lithium (solid at room temperature)"
name_Be = "Beryllium (solid at room temperature)"
name_B  = "Boron (solid at room temperature)"
name_C  = ["Carbon (amorpheous solid at room temperature)", "Carbon (graphite solid at room temperature)", "Carbon (diamond solid at room temperature)"]
name_N  = "Nitrogen  (gas N2 at room temperature)"
name_O  = "Oxygen (gas O2 at room temperature)"
name_F  = "Fluorine (halogen gas F2 at room temperature)"
name_Ne = "Neon (noble gas at room temperature)"
name_Na = "Sodium (solid at room temperature)"
name_Mg = "Magnesium (solid at room temperature)"
name_Al = "Aluminum (solid at room temperature)"
name_Si = "Silicon (solid at room temperature)"
name_P  = ["Phosphorus (white solid at room temperature)", "Phosphorus (red solid at room temperature)", "Phosphorus (violet solid at room temperature)", "Phosphorus (black solid at room temperature)"]
name_S  = ["Sulfur (solid alpha at room temperature)", "Sulfur (solid beta at room temperature)", "Sulfur (solid gamma at room temperature)"]
name_Cl = "Chlorine (halogen gas Cl2 at room temperature)"
name_Ar = "Argon (noble gas at room temperature)"
name_K  = "Potassium (solid at room temperature)"
name_Ca = "Calcium (solid at room temperature)"
name_Sc = "Scandium (solid at room temperature)"
name_Ti = "Titanium (solid at room temperature)"
name_V  = "Vanadium (solid at room temperature)"
name_Cr = "Chromium (solid at room temperature)"
name_Mn = "Manganese (solid at room temperature)"
name_Fe = "Iron (solid at room temperature)"
name_Co = "Cobalt (solid at room temperature)"
name_Ni = "Nickel (solid at room temperature)"
name_Cu = "Copper (solid at room temperature)"
name_Zn = "Zinc (solid at room temperature)"
name_Ga = "Gallium (solid at room temperature)"
name_Ge = "Germanium (solid at room temperature)"
name_As = "Arsenic (solid at room temperature)"
name_Se = ["Selenium  (gray solid at room temperature)", "Selenium  (solid alpha at room temperature)", "Selenium (vitreous solid at room temperature)"]
name_Br = "Bromine (halogen liquid Br2 at room temperature)"
name_Kr = "Krypton (noble gas at room temperature)"
name_Rb = "Rubidium (solid at room temperature)"
name_Sr = "Strontium (solid at room temperature)"
name_Y  = "Yttrium (solid at room temperature)"
name_Zr = "Zirconium (solid at room temperature)"
name_Nb = "Niobium (solid at room temperature)"
name_Mo = "Molybdenum (solid at room temperature)"
name_Tc = "technetium (solid artificially prepared at room temperature)"
name_Ru = "Ruthenium (solid at room temperature)"
name_Rh = "Rhodium (solid at room temperature)"
name_Pd = "Palladium (solid at room temperature)"
name_Ag = "Silver (solid at room temperature)"
name_Cd = "Cadmium (solid at room temperature)"
name_In = "Indium (solid at room temperature)"
name_Sn = ["Tin (gray solid alpha at room temperature)", "Tin (white solid beta at room temperature)"]
name_Sb = "Antimony (solid at room temperature)"
name_Te = "Tellurium (solid at room temperature)"
name_I  = "Iodine (halogen solid I2 at room temperature)"
name_Xe = "Xenon (noble gas at room temperature)"
name_Cs = "Cesium (solid at room temperature)"
name_Ba = "Barium (solid at room temperature)"
name_La = "Lanthanum (solid at room temperature)"
name_Ce = "Cerium (solid at room temperature)"
name_Pr = "Praseodynium (solid at room temperature)"
name_Nd = "Neodymium (solid at room temperature)"
name_Pm = "Promethium (solid artificially prepared at room temperature)"
name_Sm = "Samarium (solid at room temperature)"
name_Eu = "Europium (solid at room temperature)"
name_Gd = "Gadolinium (solid at room temperature)"
name_Tb = "Terbium (solid at room temperature)"
name_Dy = "Dysprosium (solid at room temperature)"
name_Ho = "Holmium (solid at room temperature)"
name_Er = "Erbium (solid at room temperature)"
name_Tm = "Thulium (solid at room temperature)"
name_Yb = "Ytterbium (solid at room temperature)"
name_Lu = "Lutetium (solid at room temperature)"
name_Hf = "Hafnium (solid at room temperature)"
name_Ta = "Tantalum (solid at room temperature)"
name_W  = "Tungsten (solid at room temperature)"
name_Re = "Rhenium (solid at room temperature)"
name_Os = "Osmium (solid at room temperature)"
name_Ir = "Iridium (solid at room temperature)"
name_Pt = "Platinum (solid at room temperature)"
name_Au = "Gold (solid at room temperature)"
name_Hg = "Mercury (liquid at room temperature)"
name_Tl = "Thallium (solid at room temperature)"
name_Pb = "Lead (solid at room temperature)"

##############################################
# Standard atomic weight in atomic mass unit #
##############################################
# from NIST periodic table
standard_atomic_weight_H  = 1.0079
standard_atomic_weight_He = 4.0026
standard_atomic_weight_Li = 6.94
standard_atomic_weight_Be = 9.0122
standard_atomic_weight_B  = 10.81
standard_atomic_weight_C  = 12.011
standard_atomic_weight_N  = 14.007
standard_atomic_weight_O  = 15.999
standard_atomic_weight_F  = 18.998
standard_atomic_weight_Ne = 20.180
standard_atomic_weight_Na = 22.990
standard_atomic_weight_Mg = 24.305
standard_atomic_weight_Al = 26.9815
standard_atomic_weight_Si = 28.085
standard_atomic_weight_P  = 30.974
standard_atomic_weight_S  = 32.06
standard_atomic_weight_Cl = 35.45
standard_atomic_weight_Ar = 39.948
standard_atomic_weight_K  = 39.098
standard_atomic_weight_Ca = 40.078
standard_atomic_weight_Sc = 44.956
standard_atomic_weight_Ti = 47.867
standard_atomic_weight_V  = 50.942
standard_atomic_weight_Cr = 51.996
standard_atomic_weight_Mn = 54.938
standard_atomic_weight_Fe = 55.845
standard_atomic_weight_Co = 58.993
standard_atomic_weight_Ni = 58.693
standard_atomic_weight_Cu = 63.546
standard_atomic_weight_Zn = 65.38
standard_atomic_weight_Ga = 69.723
standard_atomic_weight_Ge = 72.630
standard_atomic_weight_As = 74.922
standard_atomic_weight_Se = 78.971
standard_atomic_weight_Br = 79.904
standard_atomic_weight_Kr = 83.798
standard_atomic_weight_Rb = 85.468
standard_atomic_weight_Sr = 87.62
standard_atomic_weight_Y  = 88.906
standard_atomic_weight_Zr = 91.224
standard_atomic_weight_Nb = 92.906
standard_atomic_weight_Mo = 95.95
standard_atomic_weight_Tc = 98. 
standard_atomic_weight_Ru = 101.07
standard_atomic_weight_Rh = 102.91
standard_atomic_weight_Pd = 106.42
standard_atomic_weight_Ag = 107.87
standard_atomic_weight_Cd = 112.41
standard_atomic_weight_In = 114.82
standard_atomic_weight_Sn = 118.71
standard_atomic_weight_Sb = 121.76
standard_atomic_weight_Te = 127.60
standard_atomic_weight_I  = 126.90
standard_atomic_weight_Xe = 131.29
standard_atomic_weight_Cs = 132.91
standard_atomic_weight_Ba = 137.33
standard_atomic_weight_La = 138.91
standard_atomic_weight_Ce = 140.116
standard_atomic_weight_Pr = 140.91
standard_atomic_weight_Nd = 144.24
standard_atomic_weight_Pm = 145.
standard_atomic_weight_Sm = 150.36
standard_atomic_weight_Eu = 151.96
standard_atomic_weight_Gd = 157.25
standard_atomic_weight_Tb = 158.93
standard_atomic_weight_Dy = 162.50
standard_atomic_weight_Ho = 164.93
standard_atomic_weight_Er = 167.26
standard_atomic_weight_Tm = 168.93
standard_atomic_weight_Yb = 173.05
standard_atomic_weight_Lu = 174.97
standard_atomic_weight_Hf = 178.49
standard_atomic_weight_Ta = 180.9479
standard_atomic_weight_W  = 183.84
standard_atomic_weight_Re = 186.21
standard_atomic_weight_Os = 190.23
standard_atomic_weight_Ir = 192.22
standard_atomic_weight_Pt = 195.08
standard_atomic_weight_Au = 196.97
standard_atomic_weight_Hg = 200.59
standard_atomic_weight_Tl = 204.38
standard_atomic_weight_Pb = 207.2
#standard_atomic_weight_H2 = standard_atomic_weight_H
#standard_atomic_weight_CH = 0.5*(standard_atomic_weight_C+standard_atomic_weight_H)

######################
# Density in g/cm3 : #
######################
# from Wikipedia
density_H  = 8.3748e-5
density_He = 1.6632e-4
density_Li = 0.534
density_Be = 1.848
density_B  = 2.46
density_C  = [2.0, 1.7, 3.515]
density_N  = 1.1653e-3 # N2 
density_O  = 1.3315e-3 # O2 gas
density_F  = 1.696e-3
density_Ne = 8.3850e-4
density_Na = 0.968
density_Mg = 1.738
density_Al = 2.6989
density_Si = 2.33
density_P  = [1.823, 2.27, 2.36, 2.69] 
density_S  = [2.07, 1.96, 1.92]
density_Cl = 3.2e-3
density_Ar = 1.6620e-3
density_K  = 0.862
density_Ca = 1.55
density_Sc = 2.985
density_Ti = 4.54
density_V  = 6.0
density_Cr = 7.19
density_Mn = 7.21
density_Fe = 7.874
density_Co = 8.90
density_Ni = 8.908
density_Cu = 8.96
density_Zn = 7.14
density_Ga = 5.91
density_Ge = 5.323
density_As = 5.727
density_Se = [4.81, 4.39, 4.28]
density_Br = 3.1028
density_Kr = 3.4783e-3
density_Rb = 1.532
density_Sr = 2.64
density_Y  = 4.472
density_Zr = 6.52
density_Nb = 8.57
density_Mo = 10.22
density_Tc = 11.
density_Ru = 12.45
density_Rh = 12.41
density_Pd = 12.023
density_Ag = 10.50
density_Cd = 8.65
density_In = 7.31
density_Sn = [5.769, 7.265]
density_Sb = 6.697
density_Te = 6.24
density_I  = 4.933
density_Xe = 5.4854e-3
density_Cs = 1.93
density_Ba = 3.51
density_La = 6.162
density_Ce = 6.770
density_Pr = 6.77
density_Nd = 7.01
density_Pm = 7.26
density_Sm = 7.52
density_Eu = 5.264
density_Gd = 7.90
density_Tb = 8.23
density_Dy = 8.540
density_Ho = 8.79
density_Er = 9.066
density_Tm = 9.32
density_Yb = 6.90
density_Lu = 9.841
density_Hf = 13.31
density_Ta = 16.69
density_W  = 19.30
density_Re = 21.02
density_Os = 22.59
density_Ir = 22.56
density_Pt = 21.45
density_Au = 19.32
density_Hg = 13.534 # liquid
density_Tl = 11.85
density_Pb = 11.35

##################################
# Mean excitation energy in (eV) #
##################################
# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993)
# Tables 2.8 and 2.9
#########################################################################################################
mean_excitation_energy_H  = 19.2 # +/- 0.4
mean_excitation_energy_He = 41.8 # +/- 0.8
mean_excitation_energy_Li = 40.0 # +/- 5.
mean_excitation_energy_Be = 63.7 # +/- 3.
mean_excitation_energy_B  = 76.0 # +/- 8.
mean_excitation_energy_C  = [81.0,78.0,78.0] # [+/- 7,+/- ?,+/- ?] 
# The mean excitation energy of Diamond Carbon is assumed here to be the same as for the Graphite Carbon
mean_excitation_energy_N  = 82.0 # +/- 2.
mean_excitation_energy_O  = 95.0 # +/- 2.
mean_excitation_energy_F  = 115. # to check
mean_excitation_energy_Ne = 137. # +/- 4.
mean_excitation_energy_Na = 149. # to check
mean_excitation_energy_Mg = 156. # to check
mean_excitation_energy_Al = 166. # +/- 2.
mean_excitation_energy_Si = 173. # +/- 3.
mean_excitation_energy_P  = [173.,173.,173.,173.] # to check
# The mean excitation energy of Phosphorus is assumed here to be the same whatever its crystal structure
mean_excitation_energy_S  = [180.,180.,180.] # to check
# The mean excitation energy of Sulfur is assumed here to be the same whatever its crystal structure
mean_excitation_energy_Cl = 174. # to check
mean_excitation_energy_Ar = 188. # +/- 10.
mean_excitation_energy_K  = 190. # to check
mean_excitation_energy_Ca = 191. # +/- 8.
mean_excitation_energy_Sc = 191. # +/- 8.
mean_excitation_energy_Ti = 233. # +/- 5.
mean_excitation_energy_V  = 245. # +/- 7.
mean_excitation_energy_Cr = 257. # +/- 10.
mean_excitation_energy_Mn = 272. # +/- 10.
mean_excitation_energy_Fe = 286. # +/- 9.
mean_excitation_energy_Co = 297. # +/- 9.
mean_excitation_energy_Ni = 311. # +/- 10.
mean_excitation_energy_Cu = 322. # +/- 10.
mean_excitation_energy_Zn = 330. # +/- 10.
mean_excitation_energy_Ga = 334. # to check
mean_excitation_energy_Ge = 350. # +/- 11.
mean_excitation_energy_As = 347 # to check
mean_excitation_energy_Se = [348.,348.,348.] # to check
# The mean excitation energy of Selenium is assumed here to be the same whatever its crystal structure
mean_excitation_energy_Br = 343. # to check
mean_excitation_energy_Kr = 352. # +/- 25.
mean_excitation_energy_Rb = 363. # to check
mean_excitation_energy_Sr = 366. # to check
mean_excitation_energy_Y  = 379. # to check
mean_excitation_energy_Zr = 393. # +/- 15.
mean_excitation_energy_Nb = 417. # +/- 15.
mean_excitation_energy_Mo = 424. # +/- 15.
mean_excitation_energy_Tc = 428. # to check
mean_excitation_energy_Ru = 441. # to check
mean_excitation_energy_Rh = 449. # +/- 20.
mean_excitation_energy_Pd = 470. # +/- 20.
mean_excitation_energy_Ag = 470. # +/- 10.
mean_excitation_energy_Cd = 469. # +/- 20.
mean_excitation_energy_In = 488. # +/- 20.
mean_excitation_energy_Sn = [488.,488.] # +/- 15.
# The mean excitation energy of Tin is assumed here to be the same whatever its crystal structure
mean_excitation_energy_Sb = 487. # to check
mean_excitation_energy_Te = 485. # to check
mean_excitation_energy_I  = 474. # to check
mean_excitation_energy_Xe = 482.
mean_excitation_energy_Cs = 488. # to check
mean_excitation_energy_Ba = 491. # to check
mean_excitation_energy_La = 474. # or  501. << table 2.8
mean_excitation_energy_Ce = 508. # or  523. << table 2.8
mean_excitation_energy_Pr = 510. # or  535. << table 2.8
mean_excitation_energy_Nd = 546. # to check
mean_excitation_energy_Pm = 560. # to check
mean_excitation_energy_Sm = 561. # or  574. << table 2.8
mean_excitation_energy_Eu = 580. # to check
mean_excitation_energy_Gd = 591. # +/- 20. or 565. << table 2.9 !!! 565 << table 2.6
mean_excitation_energy_Tb = 614. # to check
mean_excitation_energy_Dy = 605. # or 628. << table 2.8
mean_excitation_energy_Ho = 640. # or 650. << table 2.8
mean_excitation_energy_Er = 658. # to check
mean_excitation_energy_Tm = 674. # to check
mean_excitation_energy_Yb = 684. # to check
mean_excitation_energy_Lu = 694. # to check
mean_excitation_energy_Hf = 671. # or 705. << table 2.8
mean_excitation_energy_Ta = 718. # +/- 30.
mean_excitation_energy_W  = 727. # +/- 30. !!!  779 << table 2.6
mean_excitation_energy_Re = 736. # to check
mean_excitation_energy_Os = 746. # to check
mean_excitation_energy_Ir = 757. # +/- 30.
mean_excitation_energy_Pt = 790. # +/- 30. !!! 786 << table 2.6
mean_excitation_energy_Au = 790. # +/- 30. !!! 790 << table 2.6
mean_excitation_energy_Hg = 800. # to check
mean_excitation_energy_Tl = 810. # to check
mean_excitation_energy_Pb = 823. # +/- 30. !!! 779 << table 2.6

################################################
#  Scaled minimum-impact parameter b used in   #
# the calculation of the the Barkas correction #
################################################
# according to Berger et al. 
# - Stopping powers and ranges for protons and alpha particles 
# -ICRU report 49 (1993) - Table 2.3

scaled_minimum_impact_parameter_H  = 0.6 #[0.6, 1.8]
scaled_minimum_impact_parameter_He = 0.6
scaled_minimum_impact_parameter_Li = 1.8
scaled_minimum_impact_parameter_Be = 1.8
scaled_minimum_impact_parameter_B  = 1.8
scaled_minimum_impact_parameter_C  = 1.8
scaled_minimum_impact_parameter_N  = 1.8
scaled_minimum_impact_parameter_O  = 1.8
scaled_minimum_impact_parameter_F  = 1.8
scaled_minimum_impact_parameter_Ne = 1.8
scaled_minimum_impact_parameter_Na = 1.4
scaled_minimum_impact_parameter_Mg = 1.4
scaled_minimum_impact_parameter_Al = 1.4
scaled_minimum_impact_parameter_Si = 1.4
scaled_minimum_impact_parameter_P  = 1.4
scaled_minimum_impact_parameter_S  = 1.4
scaled_minimum_impact_parameter_Cl = 1.4
scaled_minimum_impact_parameter_Ar = 1.8
scaled_minimum_impact_parameter_K  = 1.4
scaled_minimum_impact_parameter_Ca = 1.4
scaled_minimum_impact_parameter_Sc = 1.4
scaled_minimum_impact_parameter_Ti = 1.4
scaled_minimum_impact_parameter_V  = 1.4
scaled_minimum_impact_parameter_Cr = 1.4
scaled_minimum_impact_parameter_Mn = 1.4
scaled_minimum_impact_parameter_Fe = 1.35
scaled_minimum_impact_parameter_Co = 1.35
scaled_minimum_impact_parameter_Ni = 1.35
scaled_minimum_impact_parameter_Cu = 1.35
scaled_minimum_impact_parameter_Zn = 1.35
scaled_minimum_impact_parameter_Ga = 1.35
scaled_minimum_impact_parameter_Ge = 1.35
scaled_minimum_impact_parameter_As = 1.35
scaled_minimum_impact_parameter_Se = 1.35
scaled_minimum_impact_parameter_Br = 1.35
scaled_minimum_impact_parameter_Kr = 1.35
scaled_minimum_impact_parameter_Rb = 1.35
scaled_minimum_impact_parameter_Sr = 1.35
scaled_minimum_impact_parameter_Y  = 1.35
scaled_minimum_impact_parameter_Zr = 1.35
scaled_minimum_impact_parameter_Nb = 1.35
scaled_minimum_impact_parameter_Mo = 1.35
scaled_minimum_impact_parameter_Tc = 1.35 
scaled_minimum_impact_parameter_Ru = 1.35
scaled_minimum_impact_parameter_Rh = 1.35
scaled_minimum_impact_parameter_Pd = 1.35
scaled_minimum_impact_parameter_Ag = 1.35
scaled_minimum_impact_parameter_Cd = 1.35
scaled_minimum_impact_parameter_In = 1.35
scaled_minimum_impact_parameter_Sn = 1.35
scaled_minimum_impact_parameter_Sb = 1.3
scaled_minimum_impact_parameter_Te = 1.3
scaled_minimum_impact_parameter_I  = 1.3
scaled_minimum_impact_parameter_Xe = 1.3
scaled_minimum_impact_parameter_Cs = 1.3
scaled_minimum_impact_parameter_Ba = 1.3
scaled_minimum_impact_parameter_La = 1.3
scaled_minimum_impact_parameter_Ce = 1.3
scaled_minimum_impact_parameter_Pr = 1.3
scaled_minimum_impact_parameter_Nd = 1.3
scaled_minimum_impact_parameter_Pm = 1.3
scaled_minimum_impact_parameter_Sm = 1.3
scaled_minimum_impact_parameter_Eu = 1.3
scaled_minimum_impact_parameter_Gd = 1.3
scaled_minimum_impact_parameter_Tb = 1.3
scaled_minimum_impact_parameter_Dy = 1.3
scaled_minimum_impact_parameter_Ho = 1.3
scaled_minimum_impact_parameter_Er = 1.3
scaled_minimum_impact_parameter_Tm = 1.3
scaled_minimum_impact_parameter_Yb = 1.3
scaled_minimum_impact_parameter_Lu = 1.3
scaled_minimum_impact_parameter_Hf = 1.3
scaled_minimum_impact_parameter_Ta = 1.3
scaled_minimum_impact_parameter_W  = 1.3
scaled_minimum_impact_parameter_Re = 1.3
scaled_minimum_impact_parameter_Os = 1.3
scaled_minimum_impact_parameter_Ir = 1.3
scaled_minimum_impact_parameter_Pt = 1.3
scaled_minimum_impact_parameter_Au = 1.3
scaled_minimum_impact_parameter_Hg = 1.3
scaled_minimum_impact_parameter_Tl = 1.3
scaled_minimum_impact_parameter_Pb = 1.3

############################
# Electronic configuration #
############################

electronic_configuration_H  = "(1s)1"
electronic_configuration_He = "(1s)2"
electronic_configuration_Li = "[He]2 (2s)1"
electronic_configuration_Be = "[He]2 (2s)2"
electronic_configuration_B  = "[He]2 (2s)2 (2p)1"
electronic_configuration_C  = "[He]2 (2s)2 (2p)2"
electronic_configuration_N  = "[He]2 (2s)2 (2p)3"
electronic_configuration_O  = "[He]2 (2s)2 (2p)4"
electronic_configuration_F  = "[He]2 (2s)2 (2p)5"
electronic_configuration_Ne = "[He]2 (2s)2 (2p)6"
electronic_configuration_Na = "[Ne]10 (3s)1"
electronic_configuration_Mg = "[Ne]10 (3s)2"
electronic_configuration_Al = "[Ne]10 (3s)2 (3p)1"
electronic_configuration_Si = "[Ne]10 (3s)2 (3p)2"
electronic_configuration_P  = "[Ne]10 (3s)2 (3p)3"
electronic_configuration_S  = "[Ne]10 (3s)2 (3p)4"
electronic_configuration_Cl = "[Ne]10 (3s)2 (3p)5"
electronic_configuration_Ar = "[Ne]10 (3s)2 (3p)6"
electronic_configuration_K  = "[Ar]18 (4s)1"
electronic_configuration_Ca = "[Ar]18 (4s)2"
electronic_configuration_Sc = "[Ar]18 (4s)2 (3d)1"
electronic_configuration_Ti = "[Ar]18 (4s)2 (3d)2"
electronic_configuration_V  = "[Ar]18 (4s)2 (3d)3"
electronic_configuration_Cr = "[Ar]18 (4s)1 (3d)5"
electronic_configuration_Mn = "[Ar]18 (4s)2 (3d)5"
electronic_configuration_Fe = "[Ar]18 (4s)2 (3d)6"
electronic_configuration_Co = "[Ar]18 (4s)2 (3d)7"
electronic_configuration_Ni = "[Ar]18 (4s)2 (3d)8"
electronic_configuration_Cu = "[Ar]18 (4s)1 (3d)10"
electronic_configuration_Zn = "[Ar]18 (4s)2 (3d)10"
electronic_configuration_Ga = "[Ar]18 (4s)2 (3d)10 (4p)1"
electronic_configuration_Ge = "[Ar]18 (4s)2 (3d)10 (4p)2"
electronic_configuration_As = "[Ar]18 (4s)2 (3d)10 (4p)3"
electronic_configuration_Se = "[Ar]18 (4s)2 (3d)10 (4p)4"
electronic_configuration_Br = "[Ar]18 (4s)2 (3d)10 (4p)5"
electronic_configuration_Kr = "[Ar]18 (4s)2 (3d)10 (4p)6"
electronic_configuration_Rb = "[Kr]36 (5s)1"
electronic_configuration_Sr = "[Kr]36 (5s)2"
electronic_configuration_Y  = "[Kr]36 (5s)2 (4d)1"
electronic_configuration_Zr = "[Kr]36 (5s)2 (4d)2"
electronic_configuration_Nb = "[Kr]36 (5s)1 (4d)4"
electronic_configuration_Mo = "[Kr]36 (5s)1 (4d)5"
electronic_configuration_Tc = "[Kr]36 (5s)2 (4d)5"
electronic_configuration_Ru = "[Kr]36 (5s)1 (4d)7"
electronic_configuration_Rh = "[Kr]36 (5s)1 (4d)8"
electronic_configuration_Pd = "[Kr]36 (4d)10"
electronic_configuration_Ag = "[Kr]36 (5s)1 (4d)10"
electronic_configuration_Cd = "[Kr]36 (5s)2 (4d)10"
electronic_configuration_In = "[Kr]36 (5s)2 (4d)10 (5p)1"
electronic_configuration_Sn = "[Kr]36 (5s)2 (4d)10 (5p)2"
electronic_configuration_Sb = "[Kr]36 (5s)2 (4d)10 (5p)3"
electronic_configuration_Te = "[Kr]36 (5s)2 (4d)10 (5p)4"
electronic_configuration_I  = "[Kr]36 (5s)2 (4d)10 (5p)5"
electronic_configuration_Xe = "[Kr]36 (5s)2 (4d)10 (5p)6"
electronic_configuration_Cs = "[Xe]54 (6s)1"
electronic_configuration_Ba = "[Xe]54 (6s)2"
electronic_configuration_La = "[Xe]54 (6s)2 (5d)1"
electronic_configuration_Ce = "[Xe]54 (6s)2 (5d)1 (4f)1"
electronic_configuration_Pr = "[Xe]54 (6s)2 (4f)3"
electronic_configuration_Nd = "[Xe]54 (6s)2 (4f)4"
electronic_configuration_Pm = "[Xe]54 (6s)2 (4f)5"
electronic_configuration_Sm = "[Xe]54 (6s)2 (4f)6"
electronic_configuration_Eu = "[Xe]54 (6s)2 (4f)7"
electronic_configuration_Gd = "[Xe]54 (6s)2 (4f)7 (5d)1"
electronic_configuration_Tb = "[Xe]54 (6s)2 (4f)9"
electronic_configuration_Dy = "[Xe]54 (6s)2 (4f)10"
electronic_configuration_Ho = "[Xe]54 (6s)2 (4f)11"
electronic_configuration_Er = "[Xe]54 (6s)2 (4f)12"
electronic_configuration_Tm = "[Xe]54 (6s)2 (4f)13"
electronic_configuration_Yb = "[Xe]54 (6s)2 (4f)14"
electronic_configuration_Lu = "[Xe]54 (6s)2 (4f)14 (5d)1"
electronic_configuration_Hf = "[Xe]54 (6s)2 (4f)14 (5d)2"
electronic_configuration_Ta = "[Xe]54 (6s)2 (4f)14 (5d)3"
electronic_configuration_W  = "[Xe]54 (6s)2 (4f)14 (5d)4"
electronic_configuration_Re = "[Xe]54 (6s)2 (4f)14 (5d)5"
electronic_configuration_Os = "[Xe]54 (6s)2 (4f)14 (5d)6"
electronic_configuration_Ir = "[Xe]54 (6s)2 (4f)14 (5d)7"
electronic_configuration_Pt = "[Xe]54 (6s)1 (4f)14 (5d)9"
electronic_configuration_Au = "[Xe]54 (6s)1 (4f)14 (5d)10"
electronic_configuration_Hg = "[Xe]54 (6s)2 (4f)14 (5d)10"
electronic_configuration_Tl = "[Xe]54 (6s)2 (4f)14 (5d)10 (6p)1"
electronic_configuration_Pb = "[Xe]54 (6s)2 (4f)14 (5d)10 (6p)2"

####################
# Effective charge #
####################
# according to : _ E. Clementi and D.L. Raimondi - J. Chem. Phys. Vol. 38, No. 11 (1963)
#                _ E. Clementi, D.L. Raimondi and W.P. Reinhardt, - J. Chem. Phys. Vol. 41, No. 4 (1967)
# There are problems concerning the Aufbau exceptions La(Z=57), Ce(Z=58), Gd(Z=64), Pt(Z=78) and Au(Z=79)
# for which E. Clementi, D.L. Raimondi and W.P. Reinhardt didn't consider the correct electronic configuration
# to build the Slater's states
#                     [   1s   ,   2s   ,   2p   ,   3s   ,   3p   ,   4s   ,   3d   ,   4p   ,   5s   ,   4d   ,   5p   ,   6s   ,   4f   ,  4fbis ,   5d   ,   6p   ]
effective_charge_H  = [  1.0000]
effective_charge_He = [  1.6875]
effective_charge_Li = [  2.6906,  1.2792]
effective_charge_Be = [  3.6848,  1.9120]
effective_charge_B  = [  4.6795,  2.5762,  2.4214]
effective_charge_C  = [  5.6727,  3.2166,  3.1358]
effective_charge_N  = [  6.6651,  3.8474,  3.8340]
effective_charge_O  = [  7.6579,  4.4916,  4.4532]
effective_charge_F  = [  8.6501,  5.1276,  5.1000]
effective_charge_Ne = [  9.6421,  5.7584,  5.7584]
effective_charge_Na = [ 10.6259,  6.5714,  6.8018,  2.5074]
effective_charge_Mg = [ 11.6089,  7.3920,  7.8258,  3.3075]
effective_charge_Al = [ 12.5910,  8.2136,  8.9634,  4.1172,  4.0656]
effective_charge_Si = [ 13.5745,  9.0200,  9.9450,  4.9032,  4.2852]
effective_charge_P  = [ 14.5578,  9.8250, 10.9612,  5.6418,  4.8864]
effective_charge_S  = [ 15.5409, 10.6288, 11.9770,  6.3669,  5.4819]
effective_charge_Cl = [ 16.5239, 11.4304, 12.9932,  7.0683,  6.1161]
effective_charge_Ar = [ 17.5075, 12.2304, 14.0082,  7.7568,  6.7641]
effective_charge_K  = [ 18.4895, 13.0062, 15.0272,  8.6799,  7.7256,  3.4952]
effective_charge_Ca = [ 19.4730, 13.7764, 16.0414,  9.6015,  8.6583,  4.3980]
effective_charge_Sc = [ 20.4566, 14.5736, 17.0546, 10.3398,  9.4062,  4.6324,  7.1199]
effective_charge_Ti = [ 21.4409, 15.3766, 18.0648, 11.0331, 10.1037,  4.8168,  8.1414]
effective_charge_V  = [ 22.4256, 16.1814, 19.0728, 11.7093, 10.7850,  4.9812,  8.9829]
effective_charge_Cr = [ 23.4138, 16.9838, 20.0752, 12.3678, 11.4660,  5.1332,  9.7566]
effective_charge_Mn = [ 24.3957, 17.7938, 21.0840, 13.0179, 12.1092,  5.2832, 10.5282]
effective_charge_Fe = [ 25.3810, 18.5990, 22.0888, 13.6761, 12.7779,  5.4340, 11.1798]
effective_charge_Co = [ 26.3668, 19.4050, 23.0924, 14.3223, 13.4346,  5.5764, 11.8554]
effective_charge_Ni = [ 27.3526, 20.2126, 24.0952, 14.9610, 14.0850,  5.7108, 12.5295]
effective_charge_Cu = [ 28.3386, 21.0198, 25.0970, 15.5943, 14.7306,  5.8424, 13.2006]
effective_charge_Zn = [ 29.3245, 21.8280, 26.0980, 16.2192, 15.3693,  5.9652, 13.8783]
effective_charge_Ga = [ 30.3094, 22.5990, 27.0908, 16.9962, 16.2036,  7.0668, 15.0933, 6.2216]
effective_charge_Ge = [ 31.2937, 23.3648, 28.0822, 17.7897, 17.0136,  8.0436, 16.2513, 6.7804]
effective_charge_As = [ 32.2783, 24.1270, 29.0736, 18.5955, 17.8497,  8.9440, 17.3784, 7.4492]
effective_charge_Se = [ 33.2622, 24.8884, 30.0652, 19.4034, 18.7050,  9.7576, 18.4770, 8.2872]
effective_charge_Br = [ 34.2471, 25.6434, 31.0564, 20.2185, 19.5708, 10.5528, 19.5591, 9.0280]
effective_charge_Kr = [ 35.2316, 26.3980, 32.0470, 21.0327, 20.4342, 11.3156, 20.6259, 9.7692]
effective_charge_Rb = [ 36.2078, 27.1568, 33.0388, 21.8427, 21.3033, 12.3880, 21.6792, 10.8808,  4.9845]
effective_charge_Sr = [ 37.1911, 27.9018, 34.0304, 22.6638, 22.1676, 13.4444, 22.7262, 11.9320,  6.0705]
effective_charge_Y  = [ 38.1756, 28.6222, 35.0032, 23.5515, 23.0925, 14.2636, 25.3971, 12.7456,  6.2560, 15.9584]
effective_charge_Zr = [ 39.1590, 29.3738, 35.9928, 24.3615, 23.8455, 14.9016, 25.5669, 13.4600,  6.4455, 13.0716]
effective_charge_Nb = [ 40.1423, 30.1252, 36.9822, 25.1715, 24.6156, 15.2828, 26.2470, 14.0844,  5.9210, 11.2376]
effective_charge_Mo = [ 41.1256, 30.8768, 37.9718, 25.9815, 25.4736, 16.0964, 27.2283, 14.9768,  6.1060, 11.3924]
effective_charge_Tc = [ 42.1090, 31.6282, 38.9408, 26.7912, 26.3841, 17.1984, 28.3530, 15.8112,  7.2265, 12.8820]
effective_charge_Ru = [ 43.0923, 32.3798, 39.9508, 27.6012, 27.2211, 17.6560, 29.3589, 16.4348,  6.4845, 12.8128]
effective_charge_Rh = [ 44.0756, 33.1546, 40.9404, 28.4385, 28.1544, 18.5816, 30.4050, 17.1396,  6.6395, 13.4424]
effective_charge_Pd = [ 45.0589, 33.8828, 41.9300, 29.2212, 29.0196, 18.9860, 31.4511, 17.7232,  0.0000, 13.6176]
effective_charge_Ag = [ 46.0423, 34.6342, 42.9194, 30.0312, 29.8086, 19.8648, 32.5398, 18.5624,  6.7555, 14.7628]
effective_charge_Cd = [ 47.0256, 35.3858, 43.9090, 30.8412, 30.6915, 20.8692, 33.6069, 19.4112,  8.1920, 15.8768]
effective_charge_In = [ 48.0097, 36.1236, 44.8980, 31.6308, 31.5207, 21.7612, 34.6782, 20.3688,  9.5115, 16.9416,  8.4700]
effective_charge_Sn = [ 48.9920, 36.8594, 45.8854, 32.4198, 32.3532, 22.6580, 35.7417, 21.2652, 10.6285, 17.9700,  9.1020]
effective_charge_Sb = [ 49.9744, 37.5954, 46.8726, 33.2091, 33.1839, 23.5436, 36.7998, 22.1812, 11.6110, 18.9744,  9.9945]
effective_charge_Te = [ 50.9568, 38.3312, 47.8600, 33.9981, 34.0089, 24.4084, 37.8393, 23.1220, 12.5380, 19.9600, 10.8085]
effective_charge_I  = [ 51.9391, 39.0670, 48.8474, 34.7874, 34.8414, 25.2972, 38.9007, 24.0296, 13.4035, 20.9340, 11.6115]
effective_charge_Xe = [ 52.9215, 39.8030, 49.8346, 35.5764, 35.6676, 26.1728, 39.9468, 24.9572, 14.2180, 21.8932, 12.4245]
effective_charge_Cs = [ 53.9043, 40.5116, 50.8196, 36.3774, 36.5778, 27.0424, 40.9806, 25.8576, 15.4445, 22.8384, 13.6510,  6.3630]
effective_charge_Ba = [ 54.8861, 41.2468, 51.8096, 37.1556, 37.3164, 27.9200, 42.0243, 26.8032, 16.6195, 23.7840, 14.8005,  7.5750]
effective_charge_La = [ 55.8683, 41.9534, 52.7956, 37.9431, 38.1396, 28.7964, 43.0602, 27.7064, 17.8110, 24.7252, 15.8960,  9.3120,  1.3600, 14.0920]
effective_charge_Ce = [ 56.8481, 42.7400, 53.7824, 38.6592, 38.9595, 29.6800, 44.0853, 28.6064, 18.9135, 25.6608, 16.9655, 10.7964,  1.6760]
effective_charge_Pr = [ 57.8306, 43.4620, 54.7694, 39.5010, 39.8244, 30.3332, 45.1524, 29.0568, 17.6130, 26.2972, 15.2835,  7.7466, 21.1008,  2.0000]
effective_charge_Nd = [ 58.8132, 44.2162, 55.7566, 40.3428, 40.6890, 30.9864, 46.1568, 30.0140, 18.7425, 26.8092, 16.9610,  9.3066, 22.2660,  2.3200]
effective_charge_Pm = [ 59.7958, 44.9704, 56.7438, 41.1846, 41.5539, 31.6396, 47.0982, 30.6232, 18.8355, 27.7400, 16.4140,  9.3954, 23.1340,  2.6000]
effective_charge_Sm = [ 60.7783, 45.7348, 57.7310, 42.0264, 42.4185, 32.2924, 48.2289, 31.0880, 18.2490, 28.2396, 16.2810,  8.0118, 23.5316]
effective_charge_Eu = [ 61.7609, 46.4708, 58.7180, 42.8682, 43.2840, 32.8680, 49.2528, 31.8748, 18.5900, 28.9408, 16.5550,  8.1216, 24.3200]
effective_charge_Gd = [ 62.7435, 47.2170, 59.7054, 43.7097, 44.1492, 33.4440, 50.2770, 32.6468, 18.8820, 29.6336, 16.7640,  8.2146, 25.0136]
effective_charge_Tb = [ 63.7261, 47.9722, 60.6924, 44.5515, 45.0147, 34.0200, 51.2985, 33.3988, 19.1705, 30.3100, 16.9610,  8.3004, 25.8648]
effective_charge_Dy = [ 64.7086, 48.7094, 61.6796, 45.3933, 45.8805, 34.5920, 52.3299, 33.8260, 19.3040, 31.0180, 17.1270,  8.3436, 26.5360]
effective_charge_Ho = [ 65.6912, 49.4556, 62.6668, 46.2351, 46.7454, 35.3120, 53.3469, 34.5628, 19.5760, 31.6716, 17.3390,  8.4390, 27.4696]
effective_charge_Er = [ 66.6737, 50.2016, 63.6540, 47.0769, 47.6109, 36.2320, 54.3603, 35.1092, 19.7180, 32.2712, 17.4720,  8.4762, 27.9784]
effective_charge_Tm = [ 67.6563, 50.9478, 64.6412, 47.9184, 48.4761, 37.1376, 55.3743, 35.9880, 20.0370, 32.9440, 17.7280,  8.5842, 28.6340]
effective_charge_Yb = [ 68.6389, 51.6940, 65.6284, 48.7602, 49.3365, 37.5176, 56.3967, 36.4020, 20.1500, 33.5896, 17.8315,  8.5932, 29.4320]
effective_charge_Lu = [ 69.6195, 52.4498, 66.6110, 49.5345, 50.1663, 38.2692, 57.4188, 37.1904, 20.9550, 35.2892, 18.6800,  8.8044, 30.9312,  0.0000, 20.1130]
effective_charge_Hf = [ 70.6016, 53.1898, 67.5988, 50.3115, 50.9832, 38.9772, 58.4298, 37.9296, 21.8330, 35.5240, 19.5850,  9.1644, 32.2096,  0.0000, 16.6195]
effective_charge_Ta = [ 71.5837, 53.9298, 68.5864, 51.0915, 51.8004, 39.7588, 59.4411, 38.7348, 22.6935, 36.3240, 20.4735,  9.5250, 33.4704,  0.0000, 16.3680]
effective_charge_W  = [ 72.5657, 54.6698, 69.5742, 51.8700, 52.6176, 40.5588, 60.4524, 39.5484, 23.5415, 37.1732, 21.3255,  9.8544, 34.7108,  0.0000, 16.7420]
effective_charge_Re = [ 73.5478, 55.4098, 70.5620, 52.6485, 53.4345, 41.3564, 61.4547, 40.3732, 24.3570, 38.0544, 22.1440, 10.1160, 35.9248,  0.0000, 17.3830]
effective_charge_Os = [ 74.5299, 56.1498, 71.5498, 53.4273, 54.2517, 42.0952, 62.4747, 41.1440, 25.0950, 38.8580, 22.9100, 10.3230, 37.1528,  0.0000, 17.9970]
effective_charge_Ir = [ 75.5119, 56.8898, 72.5376, 54.2058, 55.0689, 42.8480, 63.4860, 41.9140, 25.8455, 39.7372, 23.6610, 10.5666, 38.3448,  0.0000, 18.6960]
effective_charge_Pt = [ 76.4940, 57.6298, 73.5254, 54.9843, 55.8861, 43.6388, 64.4973, 42.7304, 26.5880, 40.6300, 24.4195, 10.7514, 39.5060,  0.0000, 19.4075]
effective_charge_Au = [ 77.4761, 58.3698, 74.5132, 55.7628, 56.7030, 44.4132, 65.5083, 43.5468, 27.3275, 41.5280, 25.1700, 10.9380, 40.6496,  0.0000, 20.1265]
effective_charge_Hg = [ 78.4581, 59.1094, 75.5010, 56.5413, 57.5202, 45.2448, 66.5196, 44.4060, 28.1110, 42.4680, 25.9670, 11.1534, 41.7608,  0.0000, 20.8560]
effective_charge_Tl = [ 79.4409, 59.6842, 76.4862, 57.4191, 58.3665, 46.0788, 67.5342, 45.2168, 29.1220, 43.3888, 27.0885, 12.8196, 42.8676,  0.0000, 22.0250, 12.2538]
effective_charge_Pb = [ 80.4195, 60.4300, 77.4766, 58.1523, 59.1495, 46.8928, 68.5467, 46.0336, 30.1315, 44.3196, 28.0300, 14.1000, 43.9688,  0.0000, 23.1520, 12.3930]

######################################################
# Observed ionization potential of M-shell electrons #
######################################################
# according to J. A. Bearden and A. F. Burr, Rev. Modern Phys. Vol. 39 No.1 (1967)
# M_shell_ionization_potential_model2[i][j] = ionization potential of a MI (i=0), MII (i=1), MIII (i=2), MIV (i=3) or MV (i=4) electron
#                                             in an atom of atomic number Z = j + 1 from Hydrogen (Z=1) to Lead (Z=82)

M_shell_ionization_potential_model2 = [[   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,  17.5,  25.3,  33.9,  43.7,  53.8,  60.3,  66.5,  74.1,  83.9,  92.9, 100.7, 111.8, 119.8, 135.9, 158.1, 180.0, 203.5, 231.5, 256.5,   0. , 322.1, 357.5, 393.6, 430.3, 468.4, 504.6,   0. , 585.0, 627.1, 669.9, 717.5, 770.2, 825.6, 883.8, 943.7,1006.0,1072.1,   0. ,1217.1,1292.8,1361.3,1434.6,1511.0,1575.3,   0. ,1722.8,1800.0,1880.8,1967.5,2046.8,2128.3,2206.5,2306.8,2398.1,2491.2,2600.9,2708.0,2819.6,2931.7,3048.5,3173.7,3296.0,3424.9,3561.6,3704.1,3850.7],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   6.8,  12.4,  17.8,  25.4,  32.3,  34.6,  37.8,  42.5,  48.6,  54.0,  59.5,  68.1,  73.6,  86.6, 106.8, 127.9, 146.4, 168.2, 189.3, 222.7, 247.4, 279.8, 312.4, 344.2, 378.4, 409.7, 444.9, 482.8, 521.0, 559.1, 602.4, 650.7, 702.2, 756.4, 811.9, 869.7, 930.5, 999.0,1065.0,1136.7,1204.4,1278.8,1337.4,1402.8,1471.4,1540.7,1613.9,1688.3,1767.7,1841.8,1922.8,2005.8,2089.8,2173.0,2263.5,2365.4,2468.7,2574.9,2681.6,2792.2,2908.7,3026.5,3147.8,3278.5,3415.7,3554.2],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   6.8,  12.4,  17.8,  25.4,  32.3,  34.6,  37.8,  42.5,  48.6,  54.0,  59.5,  68.1,  73.6,  86.6, 102.9, 120.8, 140.5, 161.9, 181.5, 213.8, 238.5, 269.1, 300.3, 330.5, 363.0, 392.3, 425.0, 460.6, 496.2, 531.5, 571.4, 616.5, 664.3, 714.4, 765.6, 818.7, 874.6, 937.0, 997.6,1062.2,1123.4,1185.4,1242.2,1297.4,1356.9,1419.8,1480.6,1544.0,1611.3,1675.6,1741.2,1811.8,1884.5,1949.8,2023.6,2107.6,2194.0,2281.0,2367.3,2457.2,2550.7,2645.4,2743.0,2847.1,2956.6,3066.4],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   6.6,   3.7,   2.2,   2.3,   3.3,   3.6,   2.9,   3.6,   1.6,   8.1,  17.4,  28.7,  41.2,  56.7,  70.1,  88.9, 111.8, 135.0, 159.6, 182.4, 207.4, 230.3, 256.4, 283.6, 311.7, 340.0, 372.8, 410.5, 450.8, 493.3, 536.9, 582.5, 631.3,   0. , 739.5, 796.1, 848.5, 901.3, 951.1, 999.9,1051.5,1106.0,1160.6,1217.2,1275.0,1332.5,1391.5,1453.3,1514.6,1576.3,1639.4,1716.4,1793.2,1876.6,1948.9,2030.8,2116.1,2201.9,2291.1,2384.9,2485.1,2585.6],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   6.6,   3.7,   2.2,   2.3,   3.3,   3.6,   2.9,   3.6,   1.6,   8.1,  17.4,  28.7,  41.2,  56.7,  69.0,  88.9, 110.3, 133.1, 157.4, 180.0, 204.6, 227.0, 252.9, 279.4, 307.0, 334.7, 366.7, 403.7, 443.1, 484.8, 527.5, 572.1, 619.4, 672.3, 725.5, 780.7, 831.7, 883.3, 931.0, 977.7,1026.9,1080.2,1130.9,1185.2,1241.2,1294.9,1351.4,1409.3,1467.7,1527.8,1588.5,1661.7,1735.1,1809.2,1882.9,1960.1,2040.4,2121.6,2205.7,2294.9,2389.3,2484.0]]

######################################################
# Observed ionization potential of N-shell electrons #
######################################################
# according to J. A. Bearden and A. F. Burr, Rev. Modern Phys. Vol. 39 No.1 (1967)
# N_shell_ionization_potential_model2[i][j] = ionization potential of a NI (i=0), NII (i=1), ... or NVII (i=6) electron
#                                             in an atom of atomic number Z = j + 1 from Hydrogen (Z=1) to Lead (Z=82)

N_shell_ionization_potential_model2 = [[   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,  27.3,  24.0,  29.3,  37.7,  45.4,  51.3,  58.1,  61.8,   0. ,  74.9,  81.0,  86.4,  95.2, 107.6, 121.9, 136.5, 152.0, 168.3, 186.4,   0. , 230.8, 253.0, 270.4, 289.6, 304.5, 315.2,   0. , 345.7, 360.2, 375.8, 397.9, 416.3, 435.7, 449.1, 471.7, 487.2, 506.2, 538.1, 565.5, 595.0, 625.0, 654.3, 690.1, 722.0, 758.8, 800.3, 845.5, 893.6],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   2.5,   5.6,   5.2,  10.6,  14.8,  19.9,  25.6,  28.7,  33.9,  34.8,  38.9,  43.1,  47.9,  51.1,  62.6,  66.9,  77.4,  88.6,  98.4, 110.2, 122.7, 146.7, 172.3, 191.8, 208.8, 223.3, 263.3, 243.3, 242.0, 265.6, 283.9, 288.5, 310.2, 331.8, 343.5, 366.2, 385.9, 396.7, 410.1, 437.0, 464.8, 491.6, 517.9, 546.5, 577.1, 609.2, 643.7, 676.9, 721.3, 763.9],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   2.5,   5.6,   4.6,  10.6,  14.0,  19.9,  25.6,  28.7,  33.9,  34.8,  38.9,  43.1,  47.9,  51.1,  55.9,  66.9,  77.4,  88.6,  98.4, 110.2, 122.7, 146.7, 161.6, 179.7, 191.4, 207.2, 217.6, 224.6, 242.0, 247.4, 256.6, 270.9, 385.0, 292.9, 306.6, 320.0, 336.6, 343.5, 359.3, 380.4, 404.5, 425.3, 444.4, 468.2, 494.3, 519.0, 545.4, 571.0, 609.0, 644.5],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   3.2,   1.8,   0. ,   2.0,   2.5,   1.5,   3.3,   9.3,  16.2,  23.9,  31.4,  39.8,  49.6,   0. ,  78.8,  92.5,  98.9, 110.0, 113.2, 117.5, 120.4, 129.0, 133.2, 140.5, 147.0, 154.2, 161.0, 176.7, 179.6, 198.4, 204.8, 223.8, 241.3, 258.8, 273.7, 289.4, 311.4, 330.8, 352.0, 378.3, 406.6, 435.2],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   3.2,   1.8,   0. ,   2.0,   2.5,   1.5,   3.3,   9.3,  16.2,  23.9,  31.4,  39.8,  49.6,   0. ,  76.5,  89.9,  98.9, 110.0, 113.2, 117.5, 120.4, 129.0, 133.2, 140.5, 147.0, 154.2, 161.0, 167.6, 179.6, 184.9, 195.0, 213.7, 229.3, 245.4, 260.2, 272.8, 294.9, 313.3, 333.9, 359.8, 386.2, 412.9],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0.1,   2.0,   1.5,   0. ,   5.5,   0.0,   0.1,   2.6,   4.2,   3.7,   4.3,   5.3,   6.3,   6.9,  17.1,  25.0,  36.5,  40.6,  46.3,  63.4,  74.3,  86.4, 102.2, 122.8, 142.9],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0.1,   2.0,   1.5,   0. ,   5.5,   0.0,   0.1,   2.6,   4.2,   3.7,   4.3,   5.3,   6.3,   6.9,  17.1,  25.0,  33.6,  40.6,  46.3,  60.5,  71.1,  82.8,  98.5, 118.5, 138.1]]

######################################################
# Observed ionization potential of O-shell electrons #
######################################################
# according to J. A. Bearden and A. F. Burr, Rev. Modern Phys. Vol. 39 No.1 (1967)
# O_shell_ionization_potential_model2[i][j] = ionization potential of a OI (i=0), OII (i=1), ... or OV (i=4) electron (no electron that occupies a state OVI. OVII, OVIII or OIX > OV for Z <= 82)
#                                             in an atom of atomic number Z = j + 1 from Hydrogen (Z=1) to Lead (Z=82)

O_shell_ionization_potential_model2 = [[   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0.1,   0.9,   6.7,  11.6,  13.6,   0. ,  22.7,  39.1,  32.3,  37.8,  37.4,  37.5,   0. ,  37.4,  31.8,  36.1,  39.0,  62.9,  51.2,  59.8,  53.2,  54.1,  56.8,  64.9,  71.1,  77.1,  82.8,  83.7,  95.2, 101.7, 107.8, 120.3, 136.3, 147.3],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0.8,   1.1,   2.1,   2.3,   3.3,   0. ,  13.1,  16.6,  14.4,  19.8,  22.3,  21.1,   0. ,  21.3,  22.0,  20.3,  25.4,  26.3,  20.3,  29.4,  32.3,  23.4,  28.0,  38.1,  44.9,  46.8,  45.6,  58.0,  63.0,  65.3,  71.7,  80.5,  99.6, 104.8],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0.8,   1.1,   2.1,   2.3,   3.3,   0. ,  11.4,  14.6,  14.4,  19.8,  22.3,  21.1,   0. ,  21.3,  22.0,  20.3,  25.4,  26.3,  20.3,  29.4,  32.3,  23.4,  28.0,  30.6,  36.4,  35.6,  34.6,  45.4,  50.5,  51.7,  53.7,  57.6,  75.4,  86.0],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   5.7,   6.1,   3.5,   0. ,   3.8,   2.2,   2.5,   6.4,  15.3,  21.8],
                                       [   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   5.7,   6.1,   3.5,   0. ,   3.8,   2.2,   2.5,   6.4,  13.1,  19.2],]

##############################################################
# shell correction scale factors according to the model 1 of #
#                       Berger et al.                        #
# Stopping powers and ranges for protons and alpha particles #
#                   ICRU report 49 (1993)                    #
##############################################################

# M-shell :
M_shell_scale_parameter_model1_H  = 0.
M_shell_scale_parameter_model1_He = 0.
M_shell_scale_parameter_model1_Li = 0.
M_shell_scale_parameter_model1_Be = 0.
M_shell_scale_parameter_model1_B  = 0.
M_shell_scale_parameter_model1_C  = 0.
M_shell_scale_parameter_model1_N  = 0
M_shell_scale_parameter_model1_O  = 0.
M_shell_scale_parameter_model1_F  = 0.
M_shell_scale_parameter_model1_Ne = 0.
M_shell_scale_parameter_model1_Na = 12.0
M_shell_scale_parameter_model1_Mg = 12.0
M_shell_scale_parameter_model1_Al = 12.0
M_shell_scale_parameter_model1_Si = 12.0
M_shell_scale_parameter_model1_P  = 11.9
M_shell_scale_parameter_model1_S  = 11.7
M_shell_scale_parameter_model1_Cl = 11.5
M_shell_scale_parameter_model1_Ar = 11.2
M_shell_scale_parameter_model1_K  = 10.8
M_shell_scale_parameter_model1_Ca = 10.4
M_shell_scale_parameter_model1_Sc = 10.0
M_shell_scale_parameter_model1_Ti = 9.51
M_shell_scale_parameter_model1_V  = 8.97
M_shell_scale_parameter_model1_Cr = 8.52
M_shell_scale_parameter_model1_Mn = 8.03
M_shell_scale_parameter_model1_Fe = 7.46
M_shell_scale_parameter_model1_Co = 6.95
M_shell_scale_parameter_model1_Ni = 6.53
M_shell_scale_parameter_model1_Cu = 6.18
M_shell_scale_parameter_model1_Zn = 5.87
M_shell_scale_parameter_model1_Ga = 5.61
M_shell_scale_parameter_model1_Ge = 5.39
M_shell_scale_parameter_model1_As = 5.19
M_shell_scale_parameter_model1_Se = 5.01
M_shell_scale_parameter_model1_Br = 4.86
M_shell_scale_parameter_model1_Kr = 4.72
M_shell_scale_parameter_model1_Rb = 4.62
M_shell_scale_parameter_model1_Sr = 4.53
M_shell_scale_parameter_model1_Y  = 4.44
M_shell_scale_parameter_model1_Zr = 4.38
M_shell_scale_parameter_model1_Nb = 4.32
M_shell_scale_parameter_model1_Mo = 4.26
M_shell_scale_parameter_model1_Tc = 4.20
M_shell_scale_parameter_model1_Ru = 4.15
M_shell_scale_parameter_model1_Rh = 4.10
M_shell_scale_parameter_model1_Pd = 4.04
M_shell_scale_parameter_model1_Ag = 4.00
M_shell_scale_parameter_model1_Cd = 3.95
M_shell_scale_parameter_model1_In = 3.93
M_shell_scale_parameter_model1_Sn = 3.91
M_shell_scale_parameter_model1_Sb = 3.90
M_shell_scale_parameter_model1_Te = 3.89
M_shell_scale_parameter_model1_I  = 3.89
M_shell_scale_parameter_model1_Xe = 3.88
M_shell_scale_parameter_model1_Cs = 3.88
M_shell_scale_parameter_model1_Ba = 3.88
M_shell_scale_parameter_model1_La = 3.88
M_shell_scale_parameter_model1_Ce = 3.88
M_shell_scale_parameter_model1_Pr = 3.89
M_shell_scale_parameter_model1_Nd = 3.89
M_shell_scale_parameter_model1_Pm = 3.90
M_shell_scale_parameter_model1_Sm = 3.92
M_shell_scale_parameter_model1_Eu = 3.93

# N-shell :
N_shell_scale_parameter_model1_H  = 0.
N_shell_scale_parameter_model1_He = 0.
N_shell_scale_parameter_model1_Li = 0.
N_shell_scale_parameter_model1_Be = 0.
N_shell_scale_parameter_model1_B  = 0.
N_shell_scale_parameter_model1_C  = 0.
N_shell_scale_parameter_model1_N  = 0.
N_shell_scale_parameter_model1_O  = 0.
N_shell_scale_parameter_model1_F  = 0.
N_shell_scale_parameter_model1_Ne = 0.
N_shell_scale_parameter_model1_Na = 0.
N_shell_scale_parameter_model1_Mg = 0.
N_shell_scale_parameter_model1_Al = 0.
N_shell_scale_parameter_model1_Si = 0.
N_shell_scale_parameter_model1_P  = 0.
N_shell_scale_parameter_model1_S  = 0.
N_shell_scale_parameter_model1_Cl = 0.
N_shell_scale_parameter_model1_Ar = 0.
N_shell_scale_parameter_model1_K  = 0.
N_shell_scale_parameter_model1_Ca = 0.
N_shell_scale_parameter_model1_Sc = 0.
N_shell_scale_parameter_model1_Ti = 0.
N_shell_scale_parameter_model1_V  = 0.
N_shell_scale_parameter_model1_Cr = 0.
N_shell_scale_parameter_model1_Mn = 0.
N_shell_scale_parameter_model1_Fe = 0.
N_shell_scale_parameter_model1_Co = 0.
N_shell_scale_parameter_model1_Ni = 0.
N_shell_scale_parameter_model1_Cu = 0.
N_shell_scale_parameter_model1_Zn = 0.
N_shell_scale_parameter_model1_Ga = 0.
N_shell_scale_parameter_model1_Ge = 0.
N_shell_scale_parameter_model1_As = 75.5
N_shell_scale_parameter_model1_Se = 61.9
N_shell_scale_parameter_model1_Br = 52.2
N_shell_scale_parameter_model1_Kr = 45.1
N_shell_scale_parameter_model1_Rb = 39.6
N_shell_scale_parameter_model1_Sr = 35.4
N_shell_scale_parameter_model1_Y  = 31.9
N_shell_scale_parameter_model1_Zr = 29.1
N_shell_scale_parameter_model1_Nb = 27.2
N_shell_scale_parameter_model1_Mo = 25.8
N_shell_scale_parameter_model1_Tc = 24.5
N_shell_scale_parameter_model1_Ru = 23.6
N_shell_scale_parameter_model1_Rh = 22.7
N_shell_scale_parameter_model1_Pd = 22.0
N_shell_scale_parameter_model1_Ag = 21.4
N_shell_scale_parameter_model1_Cd = 20.9
N_shell_scale_parameter_model1_In = 20.5
N_shell_scale_parameter_model1_Sn = 20.2
N_shell_scale_parameter_model1_Sb = 19.9
N_shell_scale_parameter_model1_Te = 19.7
N_shell_scale_parameter_model1_I  = 19.5
N_shell_scale_parameter_model1_Xe = 19.3
N_shell_scale_parameter_model1_Cs = 19.2
N_shell_scale_parameter_model1_Ba = 19.1
N_shell_scale_parameter_model1_La = 18.4
N_shell_scale_parameter_model1_Ce = 18.8
N_shell_scale_parameter_model1_Pr = 18.7
N_shell_scale_parameter_model1_Nd = 18.6
N_shell_scale_parameter_model1_Pm = 18.5
N_shell_scale_parameter_model1_Sm = 18.4
N_shell_scale_parameter_model1_Eu = 18.2

############################################
# low energy empirical formula for protons #
############################################

# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993)
# Table 3.1 :

coeff_A1_low_energy_proton_empirical_formula_H  = 1.254
coeff_A1_low_energy_proton_empirical_formula_He = 1.229
coeff_A1_low_energy_proton_empirical_formula_Li = 1.411
coeff_A1_low_energy_proton_empirical_formula_Be = 2.248
coeff_A1_low_energy_proton_empirical_formula_B  = 2.474
coeff_A1_low_energy_proton_empirical_formula_C  = [2.601,2.601,2.601] # I took A2 instead of 0.000 
coeff_A1_low_energy_proton_empirical_formula_N  = 2.954
coeff_A1_low_energy_proton_empirical_formula_O  = 2.652
coeff_A1_low_energy_proton_empirical_formula_F  = 2.085
coeff_A1_low_energy_proton_empirical_formula_Ne = 1.951
coeff_A1_low_energy_proton_empirical_formula_Na = 2.542
coeff_A1_low_energy_proton_empirical_formula_Mg = 3.791
coeff_A1_low_energy_proton_empirical_formula_Al = 4.154
coeff_A1_low_energy_proton_empirical_formula_Si = 4.914
coeff_A1_low_energy_proton_empirical_formula_P  = [3.232,3.232,3.232,3.232]
coeff_A1_low_energy_proton_empirical_formula_S  = [3.447,3.447,3.447]
coeff_A1_low_energy_proton_empirical_formula_Cl = 5.301
coeff_A1_low_energy_proton_empirical_formula_Ar = 5.731
coeff_A1_low_energy_proton_empirical_formula_K  = 5.152
coeff_A1_low_energy_proton_empirical_formula_Ca = 5.521
coeff_A1_low_energy_proton_empirical_formula_Sc = 5.201
coeff_A1_low_energy_proton_empirical_formula_Ti = 4.858
coeff_A1_low_energy_proton_empirical_formula_V  = 4.479
coeff_A1_low_energy_proton_empirical_formula_Cr = 3.983
coeff_A1_low_energy_proton_empirical_formula_Mn = 3.469
coeff_A1_low_energy_proton_empirical_formula_Fe = 3.519
coeff_A1_low_energy_proton_empirical_formula_Co = 3.140
coeff_A1_low_energy_proton_empirical_formula_Ni = 3.553
coeff_A1_low_energy_proton_empirical_formula_Cu = 3.696
coeff_A1_low_energy_proton_empirical_formula_Zn = 4.210
coeff_A1_low_energy_proton_empirical_formula_Ga = 5.041
coeff_A1_low_energy_proton_empirical_formula_Ge = 5.554
coeff_A1_low_energy_proton_empirical_formula_As = 5.323
coeff_A1_low_energy_proton_empirical_formula_Se = [5.874,5.874,5.874]
coeff_A1_low_energy_proton_empirical_formula_Br = 6.658
coeff_A1_low_energy_proton_empirical_formula_Kr = 6.413
coeff_A1_low_energy_proton_empirical_formula_Rb = 5.694
coeff_A1_low_energy_proton_empirical_formula_Sr = 6.339
coeff_A1_low_energy_proton_empirical_formula_Y  = 6.407
coeff_A1_low_energy_proton_empirical_formula_Zr = 6.734
coeff_A1_low_energy_proton_empirical_formula_Nb = 6.901
coeff_A1_low_energy_proton_empirical_formula_Mo = 6.424
coeff_A1_low_energy_proton_empirical_formula_Tc = 6.799
coeff_A1_low_energy_proton_empirical_formula_Ru = 6.109
coeff_A1_low_energy_proton_empirical_formula_Rh = 5.924
coeff_A1_low_energy_proton_empirical_formula_Pd = 5.238
coeff_A1_low_energy_proton_empirical_formula_Ag = 5.345
coeff_A1_low_energy_proton_empirical_formula_Cd = 5.814
coeff_A1_low_energy_proton_empirical_formula_In = 6.229
coeff_A1_low_energy_proton_empirical_formula_Sn = [6.409,6.409]
coeff_A1_low_energy_proton_empirical_formula_Sb = 7.500
coeff_A1_low_energy_proton_empirical_formula_Te = 6.979
coeff_A1_low_energy_proton_empirical_formula_I  = 7.725
coeff_A1_low_energy_proton_empirical_formula_Xe = 8.337
coeff_A1_low_energy_proton_empirical_formula_Cs = 7.287
coeff_A1_low_energy_proton_empirical_formula_Ba = 7.899
coeff_A1_low_energy_proton_empirical_formula_La = 8.041
coeff_A1_low_energy_proton_empirical_formula_Ce = 7.488
coeff_A1_low_energy_proton_empirical_formula_Pr = 7.291
coeff_A1_low_energy_proton_empirical_formula_Nd = 7.098
coeff_A1_low_energy_proton_empirical_formula_Pm = 6.909
coeff_A1_low_energy_proton_empirical_formula_Sm = 6.728
coeff_A1_low_energy_proton_empirical_formula_Eu = 6.551
coeff_A1_low_energy_proton_empirical_formula_Gd = 6.739
coeff_A1_low_energy_proton_empirical_formula_Tb = 6.212
coeff_A1_low_energy_proton_empirical_formula_Dy = 5.517
coeff_A1_low_energy_proton_empirical_formula_Ho = 5.220
coeff_A1_low_energy_proton_empirical_formula_Er = 5.071
coeff_A1_low_energy_proton_empirical_formula_Tm = 4.926
coeff_A1_low_energy_proton_empirical_formula_Yb = 4.788
coeff_A1_low_energy_proton_empirical_formula_Lu = 4.893
coeff_A1_low_energy_proton_empirical_formula_Hf = 5.028
coeff_A1_low_energy_proton_empirical_formula_Ta = 4.738
coeff_A1_low_energy_proton_empirical_formula_W  = 4.587
coeff_A1_low_energy_proton_empirical_formula_Re = 5.201
coeff_A1_low_energy_proton_empirical_formula_Os = 5.071
coeff_A1_low_energy_proton_empirical_formula_Ir = 4.946
coeff_A1_low_energy_proton_empirical_formula_Pt = 4.477
coeff_A1_low_energy_proton_empirical_formula_Au = 4.844
coeff_A1_low_energy_proton_empirical_formula_Hg = 4.307
coeff_A1_low_energy_proton_empirical_formula_Tl = 4.723
coeff_A1_low_energy_proton_empirical_formula_Pb = 5.319

coeff_A2_low_energy_proton_empirical_formula_H  = 1.440
coeff_A2_low_energy_proton_empirical_formula_He = 1.397
coeff_A2_low_energy_proton_empirical_formula_Li = 1.600
coeff_A2_low_energy_proton_empirical_formula_Be = 2.590
coeff_A2_low_energy_proton_empirical_formula_B  = 2.815
coeff_A2_low_energy_proton_empirical_formula_C  = [2.601,2.601,2.601]
coeff_A2_low_energy_proton_empirical_formula_N  = 3.350
coeff_A2_low_energy_proton_empirical_formula_O  = 3.000
coeff_A2_low_energy_proton_empirical_formula_F  = 2.352
coeff_A2_low_energy_proton_empirical_formula_Ne = 2.199
coeff_A2_low_energy_proton_empirical_formula_Na = 2.869
coeff_A2_low_energy_proton_empirical_formula_Mg = 4.293
coeff_A2_low_energy_proton_empirical_formula_Al = 4.739
coeff_A2_low_energy_proton_empirical_formula_Si = 5.598
coeff_A2_low_energy_proton_empirical_formula_P  = [3.647,3.647,3.647,3.647]
coeff_A2_low_energy_proton_empirical_formula_S  = [3.891,3.891,3.891]
coeff_A2_low_energy_proton_empirical_formula_Cl = 6.008
coeff_A2_low_energy_proton_empirical_formula_Ar = 6.500
coeff_A2_low_energy_proton_empirical_formula_K  = 5.833
coeff_A2_low_energy_proton_empirical_formula_Ca = 6.252
coeff_A2_low_energy_proton_empirical_formula_Sc = 5.884
coeff_A2_low_energy_proton_empirical_formula_Ti = 5.489
coeff_A2_low_energy_proton_empirical_formula_V  = 5.055
coeff_A2_low_energy_proton_empirical_formula_Cr = 4.489
coeff_A2_low_energy_proton_empirical_formula_Mn = 3.907
coeff_A2_low_energy_proton_empirical_formula_Fe = 3.963
coeff_A2_low_energy_proton_empirical_formula_Co = 3.535
coeff_A2_low_energy_proton_empirical_formula_Ni = 4.004
coeff_A2_low_energy_proton_empirical_formula_Cu = 4.194
coeff_A2_low_energy_proton_empirical_formula_Zn = 4.750
coeff_A2_low_energy_proton_empirical_formula_Ga = 5.697
coeff_A2_low_energy_proton_empirical_formula_Ge = 6.300
coeff_A2_low_energy_proton_empirical_formula_As = 6.012
coeff_A2_low_energy_proton_empirical_formula_Se = [6.656,6.656,6.656]
coeff_A2_low_energy_proton_empirical_formula_Br = 7.536
coeff_A2_low_energy_proton_empirical_formula_Kr = 7.240
coeff_A2_low_energy_proton_empirical_formula_Rb = 6.429
coeff_A2_low_energy_proton_empirical_formula_Sr = 7.159
coeff_A2_low_energy_proton_empirical_formula_Y  = 7.234
coeff_A2_low_energy_proton_empirical_formula_Zr = 7.603
coeff_A2_low_energy_proton_empirical_formula_Nb = 7.791
coeff_A2_low_energy_proton_empirical_formula_Mo = 7.248
coeff_A2_low_energy_proton_empirical_formula_Tc = 7.671
coeff_A2_low_energy_proton_empirical_formula_Ru = 6.887
coeff_A2_low_energy_proton_empirical_formula_Rh = 6.677
coeff_A2_low_energy_proton_empirical_formula_Pd = 5.900
coeff_A2_low_energy_proton_empirical_formula_Ag = 6.038
coeff_A2_low_energy_proton_empirical_formula_Cd = 6.554
coeff_A2_low_energy_proton_empirical_formula_In = 7.024
coeff_A2_low_energy_proton_empirical_formula_Sn = [7.227,7.227]
coeff_A2_low_energy_proton_empirical_formula_Sb = 8.480
coeff_A2_low_energy_proton_empirical_formula_Te = 7.871
coeff_A2_low_energy_proton_empirical_formula_I  = 8.716
coeff_A2_low_energy_proton_empirical_formula_Xe = 9.425
coeff_A2_low_energy_proton_empirical_formula_Cs = 8.218
coeff_A2_low_energy_proton_empirical_formula_Ba = 8.911
coeff_A2_low_energy_proton_empirical_formula_La = 9.071
coeff_A2_low_energy_proton_empirical_formula_Ce = 8.444
coeff_A2_low_energy_proton_empirical_formula_Pr = 8.219
coeff_A2_low_energy_proton_empirical_formula_Nd = 8.000
coeff_A2_low_energy_proton_empirical_formula_Pm = 7.786
coeff_A2_low_energy_proton_empirical_formula_Sm = 7.580
coeff_A2_low_energy_proton_empirical_formula_Eu = 7.380
coeff_A2_low_energy_proton_empirical_formula_Gd = 7.592
coeff_A2_low_energy_proton_empirical_formula_Tb = 6.996
coeff_A2_low_energy_proton_empirical_formula_Dy = 6.210
coeff_A2_low_energy_proton_empirical_formula_Ho = 5.874
coeff_A2_low_energy_proton_empirical_formula_Er = 5.706
coeff_A2_low_energy_proton_empirical_formula_Tm = 5.542
coeff_A2_low_energy_proton_empirical_formula_Yb = 5.386
coeff_A2_low_energy_proton_empirical_formula_Lu = 5.505
coeff_A2_low_energy_proton_empirical_formula_Hf = 5.657
coeff_A2_low_energy_proton_empirical_formula_Ta = 5.329
coeff_A2_low_energy_proton_empirical_formula_W  = 5.160
coeff_A2_low_energy_proton_empirical_formula_Re = 5.851
coeff_A2_low_energy_proton_empirical_formula_Os = 5.704
coeff_A2_low_energy_proton_empirical_formula_Ir = 5.563
coeff_A2_low_energy_proton_empirical_formula_Pt = 5.034
coeff_A2_low_energy_proton_empirical_formula_Au = 5.458
coeff_A2_low_energy_proton_empirical_formula_Hg = 4.843
coeff_A2_low_energy_proton_empirical_formula_Tl = 5.311
coeff_A2_low_energy_proton_empirical_formula_Pb = 5.982

coeff_A3_low_energy_proton_empirical_formula_H  = 2.426e+02
coeff_A3_low_energy_proton_empirical_formula_He = 4.845e+02
coeff_A3_low_energy_proton_empirical_formula_Li = 7.256e+02
coeff_A3_low_energy_proton_empirical_formula_Be = 9.660e+02
coeff_A3_low_energy_proton_empirical_formula_B  = 1.206e+03
coeff_A3_low_energy_proton_empirical_formula_C  = [1.701e+03,1.701e+03,1.701e+03]
coeff_A3_low_energy_proton_empirical_formula_N  = 1.683e+03
coeff_A3_low_energy_proton_empirical_formula_O  = 1.920e+03
coeff_A3_low_energy_proton_empirical_formula_F  = 2.157e+03
coeff_A3_low_energy_proton_empirical_formula_Ne = 2.393e+03
coeff_A3_low_energy_proton_empirical_formula_Na = 2.628e+03
coeff_A3_low_energy_proton_empirical_formula_Mg = 2.862e+03
coeff_A3_low_energy_proton_empirical_formula_Al = 2.766e+03
coeff_A3_low_energy_proton_empirical_formula_Si = 3.193e+03
coeff_A3_low_energy_proton_empirical_formula_P  = [3.561e+03,3.561e+03,3.561e+03,3.561e+03]
coeff_A3_low_energy_proton_empirical_formula_S  = [3.792e+03,3.792e+03,3.792e+03]
coeff_A3_low_energy_proton_empirical_formula_Cl = 3.969e+03
coeff_A3_low_energy_proton_empirical_formula_Ar = 4.253e+03
coeff_A3_low_energy_proton_empirical_formula_K  = 4.482e+03
coeff_A3_low_energy_proton_empirical_formula_Ca = 4.710e+03
coeff_A3_low_energy_proton_empirical_formula_Sc = 4.938e+03
coeff_A3_low_energy_proton_empirical_formula_Ti = 5.260e+03
coeff_A3_low_energy_proton_empirical_formula_V  = 5.391e+03
coeff_A3_low_energy_proton_empirical_formula_Cr = 5.616e+03
coeff_A3_low_energy_proton_empirical_formula_Mn = 5.725e+03
coeff_A3_low_energy_proton_empirical_formula_Fe = 6.065e+03
coeff_A3_low_energy_proton_empirical_formula_Co = 6.288e+03
coeff_A3_low_energy_proton_empirical_formula_Ni = 6.205e+03
coeff_A3_low_energy_proton_empirical_formula_Cu = 4.649e+03
coeff_A3_low_energy_proton_empirical_formula_Zn = 6.953e+03
coeff_A3_low_energy_proton_empirical_formula_Ga = 7.173e+03
coeff_A3_low_energy_proton_empirical_formula_Ge = 6.496e+03
coeff_A3_low_energy_proton_empirical_formula_As = 7.611e+03
coeff_A3_low_energy_proton_empirical_formula_Se = [7.395e+03,7.395e+03,7.395e+03]
coeff_A3_low_energy_proton_empirical_formula_Br = 7.694e+03
coeff_A3_low_energy_proton_empirical_formula_Kr = 1.185e+04
coeff_A3_low_energy_proton_empirical_formula_Rb = 8.478e+03
coeff_A3_low_energy_proton_empirical_formula_Sr = 8.693e+03
coeff_A3_low_energy_proton_empirical_formula_Y  = 8.907e+03
coeff_A3_low_energy_proton_empirical_formula_Zr = 9.120e+03
coeff_A3_low_energy_proton_empirical_formula_Nb = 9.333e+03
coeff_A3_low_energy_proton_empirical_formula_Mo = 9.545e+03
coeff_A3_low_energy_proton_empirical_formula_Tc = 9.756e+03
coeff_A3_low_energy_proton_empirical_formula_Ru = 9.966e+03
coeff_A3_low_energy_proton_empirical_formula_Rh = 1.018e+04
coeff_A3_low_energy_proton_empirical_formula_Pd = 1.038e+04
coeff_A3_low_energy_proton_empirical_formula_Ag = 6.790e+03
coeff_A3_low_energy_proton_empirical_formula_Cd = 1.080e+04
coeff_A3_low_energy_proton_empirical_formula_In = 1.101e+04
coeff_A3_low_energy_proton_empirical_formula_Sn = [1.121e+04,1.121e+04]
coeff_A3_low_energy_proton_empirical_formula_Sb = 8.608e+03
coeff_A3_low_energy_proton_empirical_formula_Te = 1.162e+04
coeff_A3_low_energy_proton_empirical_formula_I  = 1.183e+04
coeff_A3_low_energy_proton_empirical_formula_Xe = 1.051e+04
coeff_A3_low_energy_proton_empirical_formula_Cs = 1.223e+04
coeff_A3_low_energy_proton_empirical_formula_Ba = 1.243e+04
coeff_A3_low_energy_proton_empirical_formula_La = 1.263e+04
coeff_A3_low_energy_proton_empirical_formula_Ce = 1.283e+04
coeff_A3_low_energy_proton_empirical_formula_Pr = 1.303e+04
coeff_A3_low_energy_proton_empirical_formula_Nd = 1.323e+04
coeff_A3_low_energy_proton_empirical_formula_Pm = 1.343e+04
coeff_A3_low_energy_proton_empirical_formula_Sm = 1.362e+04
coeff_A3_low_energy_proton_empirical_formula_Eu = 1.382e+04
coeff_A3_low_energy_proton_empirical_formula_Gd = 1.402e+04
coeff_A3_low_energy_proton_empirical_formula_Tb = 1.421e+04
coeff_A3_low_energy_proton_empirical_formula_Dy = 1.440e+04
coeff_A3_low_energy_proton_empirical_formula_Ho = 1.460e+04
coeff_A3_low_energy_proton_empirical_formula_Er = 1.479e+04
coeff_A3_low_energy_proton_empirical_formula_Tm = 1.498e+04
coeff_A3_low_energy_proton_empirical_formula_Yb = 1.517e+04
coeff_A3_low_energy_proton_empirical_formula_Lu = 1.536e+04
coeff_A3_low_energy_proton_empirical_formula_Hf = 1.555e+04
coeff_A3_low_energy_proton_empirical_formula_Ta = 1.574e+04
coeff_A3_low_energy_proton_empirical_formula_W  = 1.541e+04
coeff_A3_low_energy_proton_empirical_formula_Re = 1.612e+04
coeff_A3_low_energy_proton_empirical_formula_Os = 1.630e+04
coeff_A3_low_energy_proton_empirical_formula_Ir = 1.649e+04
coeff_A3_low_energy_proton_empirical_formula_Pt = 1.667e+04
coeff_A3_low_energy_proton_empirical_formula_Au = 7.852e+03
coeff_A3_low_energy_proton_empirical_formula_Hg = 1.704e+04
coeff_A3_low_energy_proton_empirical_formula_Tl = 1.722e+04
coeff_A3_low_energy_proton_empirical_formula_Pb = 1.740e+04

coeff_A4_low_energy_proton_empirical_formula_H  = 1.200e+04
coeff_A4_low_energy_proton_empirical_formula_He = 5.873e+03
coeff_A4_low_energy_proton_empirical_formula_Li = 3.013e+03
coeff_A4_low_energy_proton_empirical_formula_Be = 1.538e+02
coeff_A4_low_energy_proton_empirical_formula_B  = 1.060e+03
coeff_A4_low_energy_proton_empirical_formula_C  = [1.279e+03,1.279e+03,1.279e+03]
coeff_A4_low_energy_proton_empirical_formula_N  = 1.900e+03
coeff_A4_low_energy_proton_empirical_formula_O  = 2.000e+03
coeff_A4_low_energy_proton_empirical_formula_F  = 2.634e+03
coeff_A4_low_energy_proton_empirical_formula_Ne = 2.699e+03
coeff_A4_low_energy_proton_empirical_formula_Na = 1.854e+03
coeff_A4_low_energy_proton_empirical_formula_Mg = 1.009e+03
coeff_A4_low_energy_proton_empirical_formula_Al = 1.645e+02
coeff_A4_low_energy_proton_empirical_formula_Si = 2.327e+02
coeff_A4_low_energy_proton_empirical_formula_P  = [1.560e+03,1.560e+03,1.560e+03,1.560e+03]
coeff_A4_low_energy_proton_empirical_formula_S  = [1.219e+03,1.219e+03,1.219e+03]
coeff_A4_low_energy_proton_empirical_formula_Cl = 6.451e+02
coeff_A4_low_energy_proton_empirical_formula_Ar = 5.300e+02
coeff_A4_low_energy_proton_empirical_formula_K  = 5.457e+02
coeff_A4_low_energy_proton_empirical_formula_Ca = 5.533e+02
coeff_A4_low_energy_proton_empirical_formula_Sc = 5.609e+02
coeff_A4_low_energy_proton_empirical_formula_Ti = 6.511e+02
coeff_A4_low_energy_proton_empirical_formula_V  = 9.523e+02
coeff_A4_low_energy_proton_empirical_formula_Cr = 1.336e+03
coeff_A4_low_energy_proton_empirical_formula_Mn = 1.461e+03
coeff_A4_low_energy_proton_empirical_formula_Fe = 1.243e+03
coeff_A4_low_energy_proton_empirical_formula_Co = 1.372e+03
coeff_A4_low_energy_proton_empirical_formula_Ni = 5.551e+02
coeff_A4_low_energy_proton_empirical_formula_Cu = 8.113e+01
coeff_A4_low_energy_proton_empirical_formula_Zn = 2.952e+02
coeff_A4_low_energy_proton_empirical_formula_Ga = 2.026e+02
coeff_A4_low_energy_proton_empirical_formula_Ge = 1.100e+02
coeff_A4_low_energy_proton_empirical_formula_As = 2.925e+02
coeff_A4_low_energy_proton_empirical_formula_Se = [1.175e+02,1.175e+02,1.175e+02]
coeff_A4_low_energy_proton_empirical_formula_Br = 2.223e+02
coeff_A4_low_energy_proton_empirical_formula_Kr = 1.537e+02
coeff_A4_low_energy_proton_empirical_formula_Rb = 2.929e+02
coeff_A4_low_energy_proton_empirical_formula_Sr = 3.303e+02
coeff_A4_low_energy_proton_empirical_formula_Y  = 3.678e+02
coeff_A4_low_energy_proton_empirical_formula_Zr = 4.052e+02
coeff_A4_low_energy_proton_empirical_formula_Nb = 4.427e+02
coeff_A4_low_energy_proton_empirical_formula_Mo = 4.802e+02
coeff_A4_low_energy_proton_empirical_formula_Tc = 5.176e+02
coeff_A4_low_energy_proton_empirical_formula_Ru = 5.551e+02
coeff_A4_low_energy_proton_empirical_formula_Rh = 5.925e+02
coeff_A4_low_energy_proton_empirical_formula_Pd = 6.300e+02
coeff_A4_low_energy_proton_empirical_formula_Ag = 3.978e+02
coeff_A4_low_energy_proton_empirical_formula_Cd = 3.555e+02
coeff_A4_low_energy_proton_empirical_formula_In = 3.709e+02
coeff_A4_low_energy_proton_empirical_formula_Sn = [3.864e+02,3.864e+02]
coeff_A4_low_energy_proton_empirical_formula_Sb = 3.480e+02
coeff_A4_low_energy_proton_empirical_formula_Te = 3.924e+02
coeff_A4_low_energy_proton_empirical_formula_I  = 3.948e+02
coeff_A4_low_energy_proton_empirical_formula_Xe = 2.696e+02
coeff_A4_low_energy_proton_empirical_formula_Cs = 3.997e+02
coeff_A4_low_energy_proton_empirical_formula_Ba = 4.021e+02
coeff_A4_low_energy_proton_empirical_formula_La = 4.045e+02
coeff_A4_low_energy_proton_empirical_formula_Ce = 4.069e+02
coeff_A4_low_energy_proton_empirical_formula_Pr = 4.093e+02
coeff_A4_low_energy_proton_empirical_formula_Nd = 4.118e+02
coeff_A4_low_energy_proton_empirical_formula_Pm = 4.142e+02
coeff_A4_low_energy_proton_empirical_formula_Sm = 4.166e+02
coeff_A4_low_energy_proton_empirical_formula_Eu = 4.190e+02
coeff_A4_low_energy_proton_empirical_formula_Gd = 4.214e+02
coeff_A4_low_energy_proton_empirical_formula_Tb = 4.239e+02
coeff_A4_low_energy_proton_empirical_formula_Dy = 4.263e+02
coeff_A4_low_energy_proton_empirical_formula_Ho = 4.287e+02
coeff_A4_low_energy_proton_empirical_formula_Er = 4.330e+02
coeff_A4_low_energy_proton_empirical_formula_Tm = 4.335e+02
coeff_A4_low_energy_proton_empirical_formula_Yb = 4.359e+02
coeff_A4_low_energy_proton_empirical_formula_Lu = 4.384e+02
coeff_A4_low_energy_proton_empirical_formula_Hf = 4.408e+02
coeff_A4_low_energy_proton_empirical_formula_Ta = 4.432e+02
coeff_A4_low_energy_proton_empirical_formula_W  = 4.153e+02
coeff_A4_low_energy_proton_empirical_formula_Re = 4.416e+02
coeff_A4_low_energy_proton_empirical_formula_Os = 4.409e+02
coeff_A4_low_energy_proton_empirical_formula_Ir = 4.401e+02
coeff_A4_low_energy_proton_empirical_formula_Pt = 4.393e+02
coeff_A4_low_energy_proton_empirical_formula_Au = 9.758e+02
coeff_A4_low_energy_proton_empirical_formula_Hg = 4.878e+02
coeff_A4_low_energy_proton_empirical_formula_Tl = 5.370e+02
coeff_A4_low_energy_proton_empirical_formula_Pb = 5.863e+02

coeff_A5_low_energy_proton_empirical_formula_H  = 1.159e-01
coeff_A5_low_energy_proton_empirical_formula_He = 5.225e-02
coeff_A5_low_energy_proton_empirical_formula_Li = 4.578e-02
coeff_A5_low_energy_proton_empirical_formula_Be = 3.475e-02
coeff_A5_low_energy_proton_empirical_formula_B  = 2.855e-02
coeff_A5_low_energy_proton_empirical_formula_C  = [1.638e-02,1.638e-02,1.638e-02]
coeff_A5_low_energy_proton_empirical_formula_N  = 2.513e-02
coeff_A5_low_energy_proton_empirical_formula_O  = 2.230e-02
coeff_A5_low_energy_proton_empirical_formula_F  = 1.816e-02
coeff_A5_low_energy_proton_empirical_formula_Ne = 1.568e-02
coeff_A5_low_energy_proton_empirical_formula_Na = 1.472e-02
coeff_A5_low_energy_proton_empirical_formula_Mg = 1.397e-02
coeff_A5_low_energy_proton_empirical_formula_Al = 2.023e-02
coeff_A5_low_energy_proton_empirical_formula_Si = 1.419e-02
coeff_A5_low_energy_proton_empirical_formula_P  = [1.267e-02,1.267e-02,1.267e-02,1.267e-02]
coeff_A5_low_energy_proton_empirical_formula_S  = [1.211e-02,1.211e-02,1.211e-02]
coeff_A5_low_energy_proton_empirical_formula_Cl = 1.183e-02
coeff_A5_low_energy_proton_empirical_formula_Ar = 1.123e-02
coeff_A5_low_energy_proton_empirical_formula_K  = 1.129e-02
coeff_A5_low_energy_proton_empirical_formula_Ca = 1.112e-02
coeff_A5_low_energy_proton_empirical_formula_Sc = 9.995e-03
coeff_A5_low_energy_proton_empirical_formula_Ti = 8.930e-03
coeff_A5_low_energy_proton_empirical_formula_V  = 9.117e-03
coeff_A5_low_energy_proton_empirical_formula_Cr = 8.413e-03
coeff_A5_low_energy_proton_empirical_formula_Mn = 8.829e-03
coeff_A5_low_energy_proton_empirical_formula_Fe = 7.782e-03
coeff_A5_low_energy_proton_empirical_formula_Co = 7.361e-03
coeff_A5_low_energy_proton_empirical_formula_Ni = 8.763e-03
coeff_A5_low_energy_proton_empirical_formula_Cu = 2.242e-02
coeff_A5_low_energy_proton_empirical_formula_Zn = 6.809e-03
coeff_A5_low_energy_proton_empirical_formula_Ga = 6.725e-03
coeff_A5_low_energy_proton_empirical_formula_Ge = 9.689e-03
coeff_A5_low_energy_proton_empirical_formula_As = 6.447e-03
coeff_A5_low_energy_proton_empirical_formula_Se = [7.684e-03,7.684e-03,7.684e-03]
coeff_A5_low_energy_proton_empirical_formula_Br = 6.509e-03
coeff_A5_low_energy_proton_empirical_formula_Kr = 2.880e-03
coeff_A5_low_energy_proton_empirical_formula_Rb = 6.087e-03
coeff_A5_low_energy_proton_empirical_formula_Sr = 6.003e-03
coeff_A5_low_energy_proton_empirical_formula_Y  = 5.889e-03
coeff_A5_low_energy_proton_empirical_formula_Zr = 5.765e-03
coeff_A5_low_energy_proton_empirical_formula_Nb = 5.587e-03
coeff_A5_low_energy_proton_empirical_formula_Mo = 5.376e-03
coeff_A5_low_energy_proton_empirical_formula_Tc = 5.315e-03
coeff_A5_low_energy_proton_empirical_formula_Ru = 5.151e-03
coeff_A5_low_energy_proton_empirical_formula_Rh = 4.919e-03
coeff_A5_low_energy_proton_empirical_formula_Pd = 4.758e-03
coeff_A5_low_energy_proton_empirical_formula_Ag = 1.676e-02
coeff_A5_low_energy_proton_empirical_formula_Cd = 4.626e-03
coeff_A5_low_energy_proton_empirical_formula_In = 4.540e-03
coeff_A5_low_energy_proton_empirical_formula_Sn = [4.474e-03,4.474e-03]
coeff_A5_low_energy_proton_empirical_formula_Sb = 9.074e-03
coeff_A5_low_energy_proton_empirical_formula_Te = 4.402e-03
coeff_A5_low_energy_proton_empirical_formula_I  = 4.376e-03
coeff_A5_low_energy_proton_empirical_formula_Xe = 6.206e-03
coeff_A5_low_energy_proton_empirical_formula_Cs = 4.447e-03
coeff_A5_low_energy_proton_empirical_formula_Ba = 4.511e-03
coeff_A5_low_energy_proton_empirical_formula_La = 4.540e-03
coeff_A5_low_energy_proton_empirical_formula_Ce = 4.420e-03
coeff_A5_low_energy_proton_empirical_formula_Pr = 4.298e-03
coeff_A5_low_energy_proton_empirical_formula_Nd = 4.182e-03
coeff_A5_low_energy_proton_empirical_formula_Pm = 4.058e-03
coeff_A5_low_energy_proton_empirical_formula_Sm = 3.976e-03
coeff_A5_low_energy_proton_empirical_formula_Eu = 3.877e-03
coeff_A5_low_energy_proton_empirical_formula_Gd = 3.863e-03
coeff_A5_low_energy_proton_empirical_formula_Tb = 3.725e-03
coeff_A5_low_energy_proton_empirical_formula_Dy = 3.632e-03
coeff_A5_low_energy_proton_empirical_formula_Ho = 3.498e-03
coeff_A5_low_energy_proton_empirical_formula_Er = 3.405e-03
coeff_A5_low_energy_proton_empirical_formula_Tm = 3.342e-03
coeff_A5_low_energy_proton_empirical_formula_Yb = 3.292e-03
coeff_A5_low_energy_proton_empirical_formula_Lu = 3.243e-03
coeff_A5_low_energy_proton_empirical_formula_Hf = 3.195e-03
coeff_A5_low_energy_proton_empirical_formula_Ta = 3.186e-03
coeff_A5_low_energy_proton_empirical_formula_W  = 3.406e-03
coeff_A5_low_energy_proton_empirical_formula_Re = 3.122e-03
coeff_A5_low_energy_proton_empirical_formula_Os = 3.082e-03
coeff_A5_low_energy_proton_empirical_formula_Ir = 2.965e-03
coeff_A5_low_energy_proton_empirical_formula_Pt = 2.871e-03
coeff_A5_low_energy_proton_empirical_formula_Au = 2.077e-02
coeff_A5_low_energy_proton_empirical_formula_Hg = 2.882e-03
coeff_A5_low_energy_proton_empirical_formula_Tl = 2.913e-03
coeff_A5_low_energy_proton_empirical_formula_Pb = 2.871e-03

# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993) - Table 3.7 
# Rk : I imposed T1 = 3 MeV and T2 = 5 MeV for those that are not specified and
#                                           for Hydrogen   : [H]1,
#                                           for Aluminum   : [Al]13 (T1 = 8 MeV and T2 = 10 MeV are rather imposed here),
#                                           for Silicon    : [Si]14
#                                           for Titanium   : [Ti]21,
#                                           for Iron       : [Fe]26,
#                                           for Copper     : [Cu]29 (T1 = 3 MeV and T2 = 4 MeV are rather imposed here),
#                                           for Germanium  : [Ge]32,
#                                           for Krypton    : [Kr]36 (T1 = 1 MeV and T2 = 1.2 MeV are rather imposed here),
#                                           for Molybdenum : [Mo]42,
#                                           for Silver     : [Ag]47, 
#                                           for Xenon      : [Xe]54,
#                                           for Gadolinium : [Gd]64
#                                           for Tungsten   :  [W]74
#                                           for Platinum   : [Pt]78, 
#                                           for Gold       : [Au]79 (T1 = 2 MeV is rather imposed here) and
#                                           for Lead       : [Pb]82

param_T1_low_energy_proton_empirical_formula_H  = 3. #0.2
param_T1_low_energy_proton_empirical_formula_He = 0.25
param_T1_low_energy_proton_empirical_formula_Li = 3.
param_T1_low_energy_proton_empirical_formula_Be = 0.3
param_T1_low_energy_proton_empirical_formula_B  = 3.
param_T1_low_energy_proton_empirical_formula_C  = [0.2,0.2,0.2]
param_T1_low_energy_proton_empirical_formula_N  = 0.25
param_T1_low_energy_proton_empirical_formula_O  = 0.25
param_T1_low_energy_proton_empirical_formula_F  = 3.
param_T1_low_energy_proton_empirical_formula_Ne = 0.3
param_T1_low_energy_proton_empirical_formula_Na = 8.
param_T1_low_energy_proton_empirical_formula_Mg = 3.
param_T1_low_energy_proton_empirical_formula_Al = 3. #0.3
param_T1_low_energy_proton_empirical_formula_Si = 3. #0.5
param_T1_low_energy_proton_empirical_formula_P  = [3.,3.,3.,3.]
param_T1_low_energy_proton_empirical_formula_S  = [3.,3.,3.]
param_T1_low_energy_proton_empirical_formula_Cl = 3.
param_T1_low_energy_proton_empirical_formula_Ar = 0.5
param_T1_low_energy_proton_empirical_formula_K  = 3.
param_T1_low_energy_proton_empirical_formula_Ca = 3.
param_T1_low_energy_proton_empirical_formula_Sc = 3.
param_T1_low_energy_proton_empirical_formula_Ti = 3. #0.5
param_T1_low_energy_proton_empirical_formula_V  = 3.
param_T1_low_energy_proton_empirical_formula_Cr = 3.
param_T1_low_energy_proton_empirical_formula_Mn = 3.
param_T1_low_energy_proton_empirical_formula_Fe = 3. #0.5
param_T1_low_energy_proton_empirical_formula_Co = 3.
param_T1_low_energy_proton_empirical_formula_Ni = 3.
param_T1_low_energy_proton_empirical_formula_Cu = 3. #0.5
param_T1_low_energy_proton_empirical_formula_Zn = 3.
param_T1_low_energy_proton_empirical_formula_Ga = 3.
param_T1_low_energy_proton_empirical_formula_Ge = 3. #0.5
param_T1_low_energy_proton_empirical_formula_As = 3.
param_T1_low_energy_proton_empirical_formula_Se = [3.,3.,3.]
param_T1_low_energy_proton_empirical_formula_Br = 3.
param_T1_low_energy_proton_empirical_formula_Kr = 1. #0.5
param_T1_low_energy_proton_empirical_formula_Rb = 3.
param_T1_low_energy_proton_empirical_formula_Sr = 3.
param_T1_low_energy_proton_empirical_formula_Y  = 3.
param_T1_low_energy_proton_empirical_formula_Zr = 3.
param_T1_low_energy_proton_empirical_formula_Nb = 3.
param_T1_low_energy_proton_empirical_formula_Mo = 3. #0.75
param_T1_low_energy_proton_empirical_formula_Tc = 3.
param_T1_low_energy_proton_empirical_formula_Ru = 3.
param_T1_low_energy_proton_empirical_formula_Rh = 3.
param_T1_low_energy_proton_empirical_formula_Pd = 3.
param_T1_low_energy_proton_empirical_formula_Ag = 3. #0.1
param_T1_low_energy_proton_empirical_formula_Cd = 3.
param_T1_low_energy_proton_empirical_formula_In = 3.
param_T1_low_energy_proton_empirical_formula_Sn = [0.5,0.5]
param_T1_low_energy_proton_empirical_formula_Sb = 3.
param_T1_low_energy_proton_empirical_formula_Te = 3.
param_T1_low_energy_proton_empirical_formula_I  = 3.
param_T1_low_energy_proton_empirical_formula_Xe = 3. #0.5
param_T1_low_energy_proton_empirical_formula_Cs = 3.
param_T1_low_energy_proton_empirical_formula_Ba = 3.
param_T1_low_energy_proton_empirical_formula_La = 3.
param_T1_low_energy_proton_empirical_formula_Ce = 3.
param_T1_low_energy_proton_empirical_formula_Pr = 3.
param_T1_low_energy_proton_empirical_formula_Nd = 3.
param_T1_low_energy_proton_empirical_formula_Pm = 3.
param_T1_low_energy_proton_empirical_formula_Sm = 3.
param_T1_low_energy_proton_empirical_formula_Eu = 3.
param_T1_low_energy_proton_empirical_formula_Gd = 3. #0.5
param_T1_low_energy_proton_empirical_formula_Tb = 3.
param_T1_low_energy_proton_empirical_formula_Dy = 3.
param_T1_low_energy_proton_empirical_formula_Ho = 3.
param_T1_low_energy_proton_empirical_formula_Er = 3.
param_T1_low_energy_proton_empirical_formula_Tm = 3.
param_T1_low_energy_proton_empirical_formula_Yb = 3.
param_T1_low_energy_proton_empirical_formula_Lu = 3.
param_T1_low_energy_proton_empirical_formula_Hf = 3.
param_T1_low_energy_proton_empirical_formula_Ta = 3.
param_T1_low_energy_proton_empirical_formula_W  = 3. #0.3
param_T1_low_energy_proton_empirical_formula_Re = 3.
param_T1_low_energy_proton_empirical_formula_Os = 3.
param_T1_low_energy_proton_empirical_formula_Ir = 3.
param_T1_low_energy_proton_empirical_formula_Pt = 3. #0.3
param_T1_low_energy_proton_empirical_formula_Au = 2. #0.3
param_T1_low_energy_proton_empirical_formula_Hg = 3.
param_T1_low_energy_proton_empirical_formula_Tl = 3.
param_T1_low_energy_proton_empirical_formula_Pb = 3. #0.5

param_T2_low_energy_proton_empirical_formula_H  = 5. #0.5
param_T2_low_energy_proton_empirical_formula_He = 0.5
param_T2_low_energy_proton_empirical_formula_Li = 5.
param_T2_low_energy_proton_empirical_formula_Be = 0.5
param_T2_low_energy_proton_empirical_formula_B  = 5.
param_T2_low_energy_proton_empirical_formula_C  = [0.5,0.5,0.5]
param_T2_low_energy_proton_empirical_formula_N  = 0.5
param_T2_low_energy_proton_empirical_formula_O  = 0.5
param_T2_low_energy_proton_empirical_formula_F  = 5.
param_T2_low_energy_proton_empirical_formula_Ne = 1.0
param_T2_low_energy_proton_empirical_formula_Na = 10.
param_T2_low_energy_proton_empirical_formula_Mg = 5.
param_T2_low_energy_proton_empirical_formula_Al = 5. #1.0
param_T2_low_energy_proton_empirical_formula_Si = 5. #0.8
param_T2_low_energy_proton_empirical_formula_P  = [5.,5.,5.,5.]
param_T2_low_energy_proton_empirical_formula_S  = [5.,5.,5.]
param_T2_low_energy_proton_empirical_formula_Cl = 5.
param_T2_low_energy_proton_empirical_formula_Ar = 1.0
param_T2_low_energy_proton_empirical_formula_K  = 5.
param_T2_low_energy_proton_empirical_formula_Ca = 5.
param_T2_low_energy_proton_empirical_formula_Sc = 5.
param_T2_low_energy_proton_empirical_formula_Ti = 5. #1.5
param_T2_low_energy_proton_empirical_formula_V  = 5.
param_T2_low_energy_proton_empirical_formula_Cr = 5.
param_T2_low_energy_proton_empirical_formula_Mn = 5.
param_T2_low_energy_proton_empirical_formula_Fe = 5. #1.0
param_T2_low_energy_proton_empirical_formula_Co = 5.
param_T2_low_energy_proton_empirical_formula_Ni = 5.
param_T2_low_energy_proton_empirical_formula_Cu = 4. #1.0
param_T2_low_energy_proton_empirical_formula_Zn = 5.
param_T2_low_energy_proton_empirical_formula_Ga = 5.
param_T2_low_energy_proton_empirical_formula_Ge = 5. #1.0
param_T2_low_energy_proton_empirical_formula_As = 5.
param_T2_low_energy_proton_empirical_formula_Se = [5.,5.,5.]
param_T2_low_energy_proton_empirical_formula_Br = 5.
param_T2_low_energy_proton_empirical_formula_Kr = 1.2 #1.5
param_T2_low_energy_proton_empirical_formula_Rb = 5.
param_T2_low_energy_proton_empirical_formula_Sr = 5.
param_T2_low_energy_proton_empirical_formula_Y  = 5.
param_T2_low_energy_proton_empirical_formula_Zr = 5.
param_T2_low_energy_proton_empirical_formula_Nb = 5.
param_T2_low_energy_proton_empirical_formula_Mo = 5. #2.0
param_T2_low_energy_proton_empirical_formula_Tc = 5.
param_T2_low_energy_proton_empirical_formula_Ru = 5.
param_T2_low_energy_proton_empirical_formula_Rh = 5.
param_T2_low_energy_proton_empirical_formula_Pd = 5.
param_T2_low_energy_proton_empirical_formula_Ag = 5. #0.3
param_T2_low_energy_proton_empirical_formula_Cd = 5.
param_T2_low_energy_proton_empirical_formula_In = 5.
param_T2_low_energy_proton_empirical_formula_Sn = [1.5,1.5]
param_T2_low_energy_proton_empirical_formula_Sb = 5.
param_T2_low_energy_proton_empirical_formula_Te = 5.
param_T2_low_energy_proton_empirical_formula_I  = 5.
param_T2_low_energy_proton_empirical_formula_Xe = 5. #1.0
param_T2_low_energy_proton_empirical_formula_Cs = 5.
param_T2_low_energy_proton_empirical_formula_Ba = 5.
param_T2_low_energy_proton_empirical_formula_La = 5.
param_T2_low_energy_proton_empirical_formula_Ce = 5.
param_T2_low_energy_proton_empirical_formula_Pr = 5.
param_T2_low_energy_proton_empirical_formula_Nd = 5.
param_T2_low_energy_proton_empirical_formula_Pm = 5.
param_T2_low_energy_proton_empirical_formula_Sm = 5.
param_T2_low_energy_proton_empirical_formula_Eu = 5.
param_T2_low_energy_proton_empirical_formula_Gd = 5. #1.0
param_T2_low_energy_proton_empirical_formula_Tb = 5.
param_T2_low_energy_proton_empirical_formula_Dy = 5.
param_T2_low_energy_proton_empirical_formula_Ho = 5.
param_T2_low_energy_proton_empirical_formula_Er = 5.
param_T2_low_energy_proton_empirical_formula_Tm = 5.
param_T2_low_energy_proton_empirical_formula_Yb = 5.
param_T2_low_energy_proton_empirical_formula_Lu = 5.
param_T2_low_energy_proton_empirical_formula_Hf = 5.
param_T2_low_energy_proton_empirical_formula_Ta = 5.
param_T2_low_energy_proton_empirical_formula_W  = 5. #0.5
param_T2_low_energy_proton_empirical_formula_Re = 5.
param_T2_low_energy_proton_empirical_formula_Os = 5.
param_T2_low_energy_proton_empirical_formula_Ir = 5.
param_T2_low_energy_proton_empirical_formula_Pt = 5. #0.5
param_T2_low_energy_proton_empirical_formula_Au = 5. #0.5
param_T2_low_energy_proton_empirical_formula_Hg = 5.
param_T2_low_energy_proton_empirical_formula_Tl = 5.
param_T2_low_energy_proton_empirical_formula_Pb = 5. #1.0

###########################################
# low energy empirical formula for alphas #
###########################################

# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993)
# Table 3.5

coeff_a1_low_energy_alpha_empirical_formula_H  = 0.35485
coeff_a1_low_energy_alpha_empirical_formula_He = 0.58
coeff_a1_low_energy_alpha_empirical_formula_Li = 1.42
coeff_a1_low_energy_alpha_empirical_formula_Be = 2.1895
coeff_a1_low_energy_alpha_empirical_formula_B  = 3.691
coeff_a1_low_energy_alpha_empirical_formula_C  = [3.83523,3.80133,3.80133]
coeff_a1_low_energy_alpha_empirical_formula_N  = 1.9259
coeff_a1_low_energy_alpha_empirical_formula_O  = 2.81015
coeff_a1_low_energy_alpha_empirical_formula_F  = 1.533
coeff_a1_low_energy_alpha_empirical_formula_Ne = 2.303
coeff_a1_low_energy_alpha_empirical_formula_Na = 9.894
coeff_a1_low_energy_alpha_empirical_formula_Mg = 4.3
coeff_a1_low_energy_alpha_empirical_formula_Al = 2.5
coeff_a1_low_energy_alpha_empirical_formula_Si = 2.1
coeff_a1_low_energy_alpha_empirical_formula_P  = [1.729,1.729,1.729,1.729]
coeff_a1_low_energy_alpha_empirical_formula_S  = [1.402,1.402,1.402]
coeff_a1_low_energy_alpha_empirical_formula_Cl = 1.117
coeff_a1_low_energy_alpha_empirical_formula_Ar = 2.291
coeff_a1_low_energy_alpha_empirical_formula_K  = 8.554
coeff_a1_low_energy_alpha_empirical_formula_Ca = 6.297
coeff_a1_low_energy_alpha_empirical_formula_Sc = 5.307
coeff_a1_low_energy_alpha_empirical_formula_Ti = 4.71
coeff_a1_low_energy_alpha_empirical_formula_V  = 6.151
coeff_a1_low_energy_alpha_empirical_formula_Cr = 6.57
coeff_a1_low_energy_alpha_empirical_formula_Mn = 5.738
coeff_a1_low_energy_alpha_empirical_formula_Fe = 5.013
coeff_a1_low_energy_alpha_empirical_formula_Co = 4.32
coeff_a1_low_energy_alpha_empirical_formula_Ni = 4.652
coeff_a1_low_energy_alpha_empirical_formula_Cu = 3.114
coeff_a1_low_energy_alpha_empirical_formula_Zn = 3.114
coeff_a1_low_energy_alpha_empirical_formula_Ga = 3.114
coeff_a1_low_energy_alpha_empirical_formula_Ge = 5.746
coeff_a1_low_energy_alpha_empirical_formula_As = 2.792
coeff_a1_low_energy_alpha_empirical_formula_Se = [4.667,4.667,4.667]
coeff_a1_low_energy_alpha_empirical_formula_Br = 2.44
coeff_a1_low_energy_alpha_empirical_formula_Kr = 1.413
coeff_a1_low_energy_alpha_empirical_formula_Rb = 11.72
coeff_a1_low_energy_alpha_empirical_formula_Sr = 7.126
coeff_a1_low_energy_alpha_empirical_formula_Y  = 11.61
coeff_a1_low_energy_alpha_empirical_formula_Zr = 10.99
coeff_a1_low_energy_alpha_empirical_formula_Nb = 9.241
coeff_a1_low_energy_alpha_empirical_formula_Mo = 9.276
coeff_a1_low_energy_alpha_empirical_formula_Tc = 3.999
coeff_a1_low_energy_alpha_empirical_formula_Ru = 4.306
coeff_a1_low_energy_alpha_empirical_formula_Rh = 3.615
coeff_a1_low_energy_alpha_empirical_formula_Pd = 5.8
coeff_a1_low_energy_alpha_empirical_formula_Ag = 5.6
coeff_a1_low_energy_alpha_empirical_formula_Cd = 3.55
coeff_a1_low_energy_alpha_empirical_formula_In = 3.6
coeff_a1_low_energy_alpha_empirical_formula_Sn = [5.4,5.4]
coeff_a1_low_energy_alpha_empirical_formula_Sb = 3.97
coeff_a1_low_energy_alpha_empirical_formula_Te = 3.65
coeff_a1_low_energy_alpha_empirical_formula_I  = 3.118
coeff_a1_low_energy_alpha_empirical_formula_Xe = 3.949
coeff_a1_low_energy_alpha_empirical_formula_Cs = 14.4
coeff_a1_low_energy_alpha_empirical_formula_Ba = 10.99
coeff_a1_low_energy_alpha_empirical_formula_La = 16.6
coeff_a1_low_energy_alpha_empirical_formula_Ce = 10.54
coeff_a1_low_energy_alpha_empirical_formula_Pr = 10.33
coeff_a1_low_energy_alpha_empirical_formula_Nd = 10.15
coeff_a1_low_energy_alpha_empirical_formula_Pm = 9.976
coeff_a1_low_energy_alpha_empirical_formula_Sm = 9.804
coeff_a1_low_energy_alpha_empirical_formula_Eu = 14.22
coeff_a1_low_energy_alpha_empirical_formula_Gd = 9.952
coeff_a1_low_energy_alpha_empirical_formula_Tb = 9.272
coeff_a1_low_energy_alpha_empirical_formula_Dy = 10.13
coeff_a1_low_energy_alpha_empirical_formula_Ho = 8.949
coeff_a1_low_energy_alpha_empirical_formula_Er = 11.94
coeff_a1_low_energy_alpha_empirical_formula_Tm = 8.472
coeff_a1_low_energy_alpha_empirical_formula_Yb = 8.301
coeff_a1_low_energy_alpha_empirical_formula_Lu = 6.567
coeff_a1_low_energy_alpha_empirical_formula_Hf = 5.951
coeff_a1_low_energy_alpha_empirical_formula_Ta = 7.495
coeff_a1_low_energy_alpha_empirical_formula_W  = 6.335
coeff_a1_low_energy_alpha_empirical_formula_Re = 4.314
coeff_a1_low_energy_alpha_empirical_formula_Os = 4.02
coeff_a1_low_energy_alpha_empirical_formula_Ir = 3.836
coeff_a1_low_energy_alpha_empirical_formula_Pt = 4.68
coeff_a1_low_energy_alpha_empirical_formula_Au = 3.223
coeff_a1_low_energy_alpha_empirical_formula_Hg = 2.892
coeff_a1_low_energy_alpha_empirical_formula_Tl = 4.728
coeff_a1_low_energy_alpha_empirical_formula_Pb = 6.18

coeff_a2_low_energy_alpha_empirical_formula_H  = 0.6456
coeff_a2_low_energy_alpha_empirical_formula_He = 0.59
coeff_a2_low_energy_alpha_empirical_formula_Li = 0.49
coeff_a2_low_energy_alpha_empirical_formula_Be = 0.47183
coeff_a2_low_energy_alpha_empirical_formula_B  = 0.4128
coeff_a2_low_energy_alpha_empirical_formula_C  = [0.42993,0.41590,0.41590]
coeff_a2_low_energy_alpha_empirical_formula_N  = 0.5550
coeff_a2_low_energy_alpha_empirical_formula_O  = 0.4759
coeff_a2_low_energy_alpha_empirical_formula_F  = 0.531
coeff_a2_low_energy_alpha_empirical_formula_Ne = 0.4861
coeff_a2_low_energy_alpha_empirical_formula_Na = 0.3081
coeff_a2_low_energy_alpha_empirical_formula_Mg = 0.47
coeff_a2_low_energy_alpha_empirical_formula_Al = 0.625
coeff_a2_low_energy_alpha_empirical_formula_Si = 0.65
coeff_a2_low_energy_alpha_empirical_formula_P  = [0.6562,0.6562,0.6562,0.6562]
coeff_a2_low_energy_alpha_empirical_formula_S  = [0.6791,0.6791,0.6791]
coeff_a2_low_energy_alpha_empirical_formula_Cl = 0.7044
coeff_a2_low_energy_alpha_empirical_formula_Ar = 0.6284
coeff_a2_low_energy_alpha_empirical_formula_K  = 0.3817
coeff_a2_low_energy_alpha_empirical_formula_Ca = 0.4622
coeff_a2_low_energy_alpha_empirical_formula_Sc = 0.4918
coeff_a2_low_energy_alpha_empirical_formula_Ti = 0.5087
coeff_a2_low_energy_alpha_empirical_formula_V  = 0.4524
coeff_a2_low_energy_alpha_empirical_formula_Cr = 0.4322
coeff_a2_low_energy_alpha_empirical_formula_Mn = 0.4492
coeff_a2_low_energy_alpha_empirical_formula_Fe = 0.4707
coeff_a2_low_energy_alpha_empirical_formula_Co = 0.4947
coeff_a2_low_energy_alpha_empirical_formula_Ni = 0.4571
coeff_a2_low_energy_alpha_empirical_formula_Cu = 0.5236
coeff_a2_low_energy_alpha_empirical_formula_Zn = 0.5236
coeff_a2_low_energy_alpha_empirical_formula_Ga = 0.5236
coeff_a2_low_energy_alpha_empirical_formula_Ge = 0.4662
coeff_a2_low_energy_alpha_empirical_formula_As = 0.6346
coeff_a2_low_energy_alpha_empirical_formula_Se = [0.5095,0.5095,0.5095]
coeff_a2_low_energy_alpha_empirical_formula_Br = 0.6346
coeff_a2_low_energy_alpha_empirical_formula_Kr = 0.7377
coeff_a2_low_energy_alpha_empirical_formula_Rb = 0.3826
coeff_a2_low_energy_alpha_empirical_formula_Sr = 0.4804
coeff_a2_low_energy_alpha_empirical_formula_Y  = 0.3955
coeff_a2_low_energy_alpha_empirical_formula_Zr = 0.41
coeff_a2_low_energy_alpha_empirical_formula_Nb = 0.4275
coeff_a2_low_energy_alpha_empirical_formula_Mo = 0.418
coeff_a2_low_energy_alpha_empirical_formula_Tc = 0.6152
coeff_a2_low_energy_alpha_empirical_formula_Ru = 0.5658
coeff_a2_low_energy_alpha_empirical_formula_Rh = 0.6197
coeff_a2_low_energy_alpha_empirical_formula_Pd = 0.49
coeff_a2_low_energy_alpha_empirical_formula_Ag = 0.49
coeff_a2_low_energy_alpha_empirical_formula_Cd = 0.6068
coeff_a2_low_energy_alpha_empirical_formula_In = 0.62
coeff_a2_low_energy_alpha_empirical_formula_Sn = [0.53,0.53]
coeff_a2_low_energy_alpha_empirical_formula_Sb = 0.6459
coeff_a2_low_energy_alpha_empirical_formula_Te = 0.64
coeff_a2_low_energy_alpha_empirical_formula_I  = 0.6519
coeff_a2_low_energy_alpha_empirical_formula_Xe = 0.6209
coeff_a2_low_energy_alpha_empirical_formula_Cs = 0.3923
coeff_a2_low_energy_alpha_empirical_formula_Ba = 0.4599
coeff_a2_low_energy_alpha_empirical_formula_La = 0.3773
coeff_a2_low_energy_alpha_empirical_formula_Ce = 0.4533
coeff_a2_low_energy_alpha_empirical_formula_Pr = 0.4502
coeff_a2_low_energy_alpha_empirical_formula_Nd = 0.4471
coeff_a2_low_energy_alpha_empirical_formula_Pm = 0.4439
coeff_a2_low_energy_alpha_empirical_formula_Sm = 0.4408
coeff_a2_low_energy_alpha_empirical_formula_Eu = 0.363
coeff_a2_low_energy_alpha_empirical_formula_Gd = 0.4318
coeff_a2_low_energy_alpha_empirical_formula_Tb = 0.4345
coeff_a2_low_energy_alpha_empirical_formula_Dy = 0.4146
coeff_a2_low_energy_alpha_empirical_formula_Ho = 0.4304
coeff_a2_low_energy_alpha_empirical_formula_Er = 0.3783
coeff_a2_low_energy_alpha_empirical_formula_Tm = 0.4405
coeff_a2_low_energy_alpha_empirical_formula_Yb = 0.4399
coeff_a2_low_energy_alpha_empirical_formula_Lu = 0.4858
coeff_a2_low_energy_alpha_empirical_formula_Hf = 0.5016
coeff_a2_low_energy_alpha_empirical_formula_Ta = 0.4523
coeff_a2_low_energy_alpha_empirical_formula_W  = 0.4825
coeff_a2_low_energy_alpha_empirical_formula_Re = 0.5558
coeff_a2_low_energy_alpha_empirical_formula_Os = 0.5681
coeff_a2_low_energy_alpha_empirical_formula_Ir = 0.5765
coeff_a2_low_energy_alpha_empirical_formula_Pt = 0.5247
coeff_a2_low_energy_alpha_empirical_formula_Au = 0.5883
coeff_a2_low_energy_alpha_empirical_formula_Hg = 0.6204
coeff_a2_low_energy_alpha_empirical_formula_Tl = 0.5522
coeff_a2_low_energy_alpha_empirical_formula_Pb = 0.52

coeff_a3_low_energy_alpha_empirical_formula_H  = 6.01525
coeff_a3_low_energy_alpha_empirical_formula_He = 6.3
coeff_a3_low_energy_alpha_empirical_formula_Li = 12.25
coeff_a3_low_energy_alpha_empirical_formula_Be = 7.2362
coeff_a3_low_energy_alpha_empirical_formula_B  = 18.48
coeff_a3_low_energy_alpha_empirical_formula_C  = [12.6125,12.9966,12.9966]
coeff_a3_low_energy_alpha_empirical_formula_N  = 27.15125
coeff_a3_low_energy_alpha_empirical_formula_O  = 50.0253
coeff_a3_low_energy_alpha_empirical_formula_F  = 40.44
coeff_a3_low_energy_alpha_empirical_formula_Ne = 37.01
coeff_a3_low_energy_alpha_empirical_formula_Na = 23.65
coeff_a3_low_energy_alpha_empirical_formula_Mg = 34.3
coeff_a3_low_energy_alpha_empirical_formula_Al = 45.7
coeff_a3_low_energy_alpha_empirical_formula_Si = 49.34
coeff_a3_low_energy_alpha_empirical_formula_P  = [53.41,53.41,53.41,53.41]
coeff_a3_low_energy_alpha_empirical_formula_S  = [58.98,58.98,58.98]
coeff_a3_low_energy_alpha_empirical_formula_Cl = 69.69
coeff_a3_low_energy_alpha_empirical_formula_Ar = 73.88
coeff_a3_low_energy_alpha_empirical_formula_K  = 83.61
coeff_a3_low_energy_alpha_empirical_formula_Ca = 65.39
coeff_a3_low_energy_alpha_empirical_formula_Sc = 61.74
coeff_a3_low_energy_alpha_empirical_formula_Ti = 65.28
coeff_a3_low_energy_alpha_empirical_formula_V  = 83.0
coeff_a3_low_energy_alpha_empirical_formula_Cr = 84.76
coeff_a3_low_energy_alpha_empirical_formula_Mn = 84.6
coeff_a3_low_energy_alpha_empirical_formula_Fe = 85.8
coeff_a3_low_energy_alpha_empirical_formula_Co = 76.14
coeff_a3_low_energy_alpha_empirical_formula_Ni = 80.73
coeff_a3_low_energy_alpha_empirical_formula_Cu = 76.67
coeff_a3_low_energy_alpha_empirical_formula_Zn = 76.67
coeff_a3_low_energy_alpha_empirical_formula_Ga = 76.67
coeff_a3_low_energy_alpha_empirical_formula_Ge = 79.24
coeff_a3_low_energy_alpha_empirical_formula_As = 106.1
coeff_a3_low_energy_alpha_empirical_formula_Se = [124.3,124.3,124.3]
coeff_a3_low_energy_alpha_empirical_formula_Br = 105.0
coeff_a3_low_energy_alpha_empirical_formula_Kr = 147.9
coeff_a3_low_energy_alpha_empirical_formula_Rb = 102.8
coeff_a3_low_energy_alpha_empirical_formula_Sr = 119.3
coeff_a3_low_energy_alpha_empirical_formula_Y  = 146.7
coeff_a3_low_energy_alpha_empirical_formula_Zr = 163.9
coeff_a3_low_energy_alpha_empirical_formula_Nb = 163.1
coeff_a3_low_energy_alpha_empirical_formula_Mo = 157.1
coeff_a3_low_energy_alpha_empirical_formula_Tc = 97.6
coeff_a3_low_energy_alpha_empirical_formula_Ru = 97.99
coeff_a3_low_energy_alpha_empirical_formula_Rh = 86.26
coeff_a3_low_energy_alpha_empirical_formula_Pd = 147.2
coeff_a3_low_energy_alpha_empirical_formula_Ag = 130.0
coeff_a3_low_energy_alpha_empirical_formula_Cd = 124.7
coeff_a3_low_energy_alpha_empirical_formula_In = 105.8
coeff_a3_low_energy_alpha_empirical_formula_Sn = [103.1,103.1]
coeff_a3_low_energy_alpha_empirical_formula_Sb = 131.8
coeff_a3_low_energy_alpha_empirical_formula_Te = 126.8
coeff_a3_low_energy_alpha_empirical_formula_I  = 164.9
coeff_a3_low_energy_alpha_empirical_formula_Xe = 200.5
coeff_a3_low_energy_alpha_empirical_formula_Cs = 152.5
coeff_a3_low_energy_alpha_empirical_formula_Ba = 138.4
coeff_a3_low_energy_alpha_empirical_formula_La = 224.1
coeff_a3_low_energy_alpha_empirical_formula_Ce = 159.3
coeff_a3_low_energy_alpha_empirical_formula_Pr = 162.0
coeff_a3_low_energy_alpha_empirical_formula_Nd = 165.6
coeff_a3_low_energy_alpha_empirical_formula_Pm = 168.0
coeff_a3_low_energy_alpha_empirical_formula_Sm = 176.2
coeff_a3_low_energy_alpha_empirical_formula_Eu = 228.4
coeff_a3_low_energy_alpha_empirical_formula_Gd = 233.5
coeff_a3_low_energy_alpha_empirical_formula_Tb = 210.0
coeff_a3_low_energy_alpha_empirical_formula_Dy = 225.7
coeff_a3_low_energy_alpha_empirical_formula_Ho = 213.3
coeff_a3_low_energy_alpha_empirical_formula_Er = 247.2
coeff_a3_low_energy_alpha_empirical_formula_Tm = 195.5
coeff_a3_low_energy_alpha_empirical_formula_Yb = 203.7
coeff_a3_low_energy_alpha_empirical_formula_Lu = 193.0
coeff_a3_low_energy_alpha_empirical_formula_Hf = 196.1
coeff_a3_low_energy_alpha_empirical_formula_Ta = 251.4
coeff_a3_low_energy_alpha_empirical_formula_W  = 255.1
coeff_a3_low_energy_alpha_empirical_formula_Re = 214.8
coeff_a3_low_energy_alpha_empirical_formula_Os = 219.9
coeff_a3_low_energy_alpha_empirical_formula_Ir = 210.2
coeff_a3_low_energy_alpha_empirical_formula_Pt = 244.7
coeff_a3_low_energy_alpha_empirical_formula_Au = 232.7
coeff_a3_low_energy_alpha_empirical_formula_Hg = 208.6
coeff_a3_low_energy_alpha_empirical_formula_Tl = 217.0
coeff_a3_low_energy_alpha_empirical_formula_Pb = 170.0

coeff_a4_low_energy_alpha_empirical_formula_H  = 20.8933
coeff_a4_low_energy_alpha_empirical_formula_He = 130.0
coeff_a4_low_energy_alpha_empirical_formula_Li = 32.0
coeff_a4_low_energy_alpha_empirical_formula_Be = 134.30
coeff_a4_low_energy_alpha_empirical_formula_B  = 50.72
coeff_a4_low_energy_alpha_empirical_formula_C  = [227.41,117.83,117.83]
coeff_a4_low_energy_alpha_empirical_formula_N  = 26.0665
coeff_a4_low_energy_alpha_empirical_formula_O  = 10.556
coeff_a4_low_energy_alpha_empirical_formula_F  = 18.41
coeff_a4_low_energy_alpha_empirical_formula_Ne = 37.96
coeff_a4_low_energy_alpha_empirical_formula_Na = 0.384
coeff_a4_low_energy_alpha_empirical_formula_Mg = 3.3
coeff_a4_low_energy_alpha_empirical_formula_Al = 0.1
coeff_a4_low_energy_alpha_empirical_formula_Si = 1.788
coeff_a4_low_energy_alpha_empirical_formula_P  = [2.405,2.405,2.405,2.405]
coeff_a4_low_energy_alpha_empirical_formula_S  = [3.528,3.528,3.528]
coeff_a4_low_energy_alpha_empirical_formula_Cl = 3.705
coeff_a4_low_energy_alpha_empirical_formula_Ar = 4.478
coeff_a4_low_energy_alpha_empirical_formula_K  = 11.84
coeff_a4_low_energy_alpha_empirical_formula_Ca = 10.14
coeff_a4_low_energy_alpha_empirical_formula_Sc = 12.4
coeff_a4_low_energy_alpha_empirical_formula_Ti = 8.806
coeff_a4_low_energy_alpha_empirical_formula_V  = 18.31
coeff_a4_low_energy_alpha_empirical_formula_Cr = 15.53
coeff_a4_low_energy_alpha_empirical_formula_Mn = 14.18
coeff_a4_low_energy_alpha_empirical_formula_Fe = 16.55
coeff_a4_low_energy_alpha_empirical_formula_Co = 10.85
coeff_a4_low_energy_alpha_empirical_formula_Ni = 22.0
coeff_a4_low_energy_alpha_empirical_formula_Cu = 7.62
coeff_a4_low_energy_alpha_empirical_formula_Zn = 7.62
coeff_a4_low_energy_alpha_empirical_formula_Ga = 7.62
coeff_a4_low_energy_alpha_empirical_formula_Ge = 1.185
coeff_a4_low_energy_alpha_empirical_formula_As = 0.2986
coeff_a4_low_energy_alpha_empirical_formula_Se = [2.102,2.102,2.102]
coeff_a4_low_energy_alpha_empirical_formula_Br = 0.83
coeff_a4_low_energy_alpha_empirical_formula_Kr = 1.466
coeff_a4_low_energy_alpha_empirical_formula_Rb = 9.231
coeff_a4_low_energy_alpha_empirical_formula_Sr = 5.784
coeff_a4_low_energy_alpha_empirical_formula_Y  = 7.031
coeff_a4_low_energy_alpha_empirical_formula_Zr = 7.1
coeff_a4_low_energy_alpha_empirical_formula_Nb = 7.954
coeff_a4_low_energy_alpha_empirical_formula_Mo = 8.038
coeff_a4_low_energy_alpha_empirical_formula_Tc = 1.297
coeff_a4_low_energy_alpha_empirical_formula_Ru = 5.514
coeff_a4_low_energy_alpha_empirical_formula_Rh = 0.333
coeff_a4_low_energy_alpha_empirical_formula_Pd = 6.903
coeff_a4_low_energy_alpha_empirical_formula_Ag = 10.0
coeff_a4_low_energy_alpha_empirical_formula_Cd = 1.112
coeff_a4_low_energy_alpha_empirical_formula_In = 0.1692
coeff_a4_low_energy_alpha_empirical_formula_Sn = [3.931,3.931]
coeff_a4_low_energy_alpha_empirical_formula_Sb = 0.2233
coeff_a4_low_energy_alpha_empirical_formula_Te = 0.6834
coeff_a4_low_energy_alpha_empirical_formula_I  = 1.208
coeff_a4_low_energy_alpha_empirical_formula_Xe = 1.878
coeff_a4_low_energy_alpha_empirical_formula_Cs = 8.354
coeff_a4_low_energy_alpha_empirical_formula_Ba = 4.811
coeff_a4_low_energy_alpha_empirical_formula_La = 6.28
coeff_a4_low_energy_alpha_empirical_formula_Ce = 4.832
coeff_a4_low_energy_alpha_empirical_formula_Pr = 5.132
coeff_a4_low_energy_alpha_empirical_formula_Nd = 5.378
coeff_a4_low_energy_alpha_empirical_formula_Pm = 5.721
coeff_a4_low_energy_alpha_empirical_formula_Sm = 5.675
coeff_a4_low_energy_alpha_empirical_formula_Eu = 7.024
coeff_a4_low_energy_alpha_empirical_formula_Gd = 5.065
coeff_a4_low_energy_alpha_empirical_formula_Tb = 4.911
coeff_a4_low_energy_alpha_empirical_formula_Dy = 5.525
coeff_a4_low_energy_alpha_empirical_formula_Ho = 5.071
coeff_a4_low_energy_alpha_empirical_formula_Er = 6.655
coeff_a4_low_energy_alpha_empirical_formula_Tm = 4.051
coeff_a4_low_energy_alpha_empirical_formula_Yb = 3.667
coeff_a4_low_energy_alpha_empirical_formula_Lu = 2.65
coeff_a4_low_energy_alpha_empirical_formula_Hf = 2.662
coeff_a4_low_energy_alpha_empirical_formula_Ta = 3.433
coeff_a4_low_energy_alpha_empirical_formula_W  = 2.834
coeff_a4_low_energy_alpha_empirical_formula_Re = 2.354
coeff_a4_low_energy_alpha_empirical_formula_Os = 2.402
coeff_a4_low_energy_alpha_empirical_formula_Ir = 2.742
coeff_a4_low_energy_alpha_empirical_formula_Pt = 2.749
coeff_a4_low_energy_alpha_empirical_formula_Au = 2.954
coeff_a4_low_energy_alpha_empirical_formula_Hg = 2.415
coeff_a4_low_energy_alpha_empirical_formula_Tl = 3.091
coeff_a4_low_energy_alpha_empirical_formula_Pb = 4.0

coeff_a5_low_energy_alpha_empirical_formula_H  = 4.3515
coeff_a5_low_energy_alpha_empirical_formula_He = 44.07
coeff_a5_low_energy_alpha_empirical_formula_Li = 9.161
coeff_a5_low_energy_alpha_empirical_formula_Be = 197.96
coeff_a5_low_energy_alpha_empirical_formula_B  = 9.0
coeff_a5_low_energy_alpha_empirical_formula_C  = [188.97,242.28,242.28]
coeff_a5_low_energy_alpha_empirical_formula_N  = 6.2768
coeff_a5_low_energy_alpha_empirical_formula_O  = 1.0382
coeff_a5_low_energy_alpha_empirical_formula_F  = 2.718
coeff_a5_low_energy_alpha_empirical_formula_Ne = 5.092
coeff_a5_low_energy_alpha_empirical_formula_Na = 92.93
coeff_a5_low_energy_alpha_empirical_formula_Mg = 12.74
coeff_a5_low_energy_alpha_empirical_formula_Al = 4.359
coeff_a5_low_energy_alpha_empirical_formula_Si = 4.133
coeff_a5_low_energy_alpha_empirical_formula_P  = [3.845,3.845,3.845,3.845]
coeff_a5_low_energy_alpha_empirical_formula_S  = [3.211,3.211,3.211]
coeff_a5_low_energy_alpha_empirical_formula_Cl = 2.156
coeff_a5_low_energy_alpha_empirical_formula_Ar = 2.066
coeff_a5_low_energy_alpha_empirical_formula_K  = 1.875
coeff_a5_low_energy_alpha_empirical_formula_Ca = 5.036
coeff_a5_low_energy_alpha_empirical_formula_Sc = 6.665
coeff_a5_low_energy_alpha_empirical_formula_Ti = 5.948
coeff_a5_low_energy_alpha_empirical_formula_V  = 2.71
coeff_a5_low_energy_alpha_empirical_formula_Cr = 2.779
coeff_a5_low_energy_alpha_empirical_formula_Mn = 3.101
coeff_a5_low_energy_alpha_empirical_formula_Fe = 3.211
coeff_a5_low_energy_alpha_empirical_formula_Co = 5.441
coeff_a5_low_energy_alpha_empirical_formula_Ni = 4.952
coeff_a5_low_energy_alpha_empirical_formula_Cu = 6.385
coeff_a5_low_energy_alpha_empirical_formula_Zn = 7.502
coeff_a5_low_energy_alpha_empirical_formula_Ga = 8.514
coeff_a5_low_energy_alpha_empirical_formula_Ge = 7.993
coeff_a5_low_energy_alpha_empirical_formula_As = 2.331
coeff_a5_low_energy_alpha_empirical_formula_Se = [1.667,1.667,1.667]
coeff_a5_low_energy_alpha_empirical_formula_Br = 2.851
coeff_a5_low_energy_alpha_empirical_formula_Kr = 1.016
coeff_a5_low_energy_alpha_empirical_formula_Rb = 4.371
coeff_a5_low_energy_alpha_empirical_formula_Sr = 2.454
coeff_a5_low_energy_alpha_empirical_formula_Y  = 1.423
coeff_a5_low_energy_alpha_empirical_formula_Zr = 1.052
coeff_a5_low_energy_alpha_empirical_formula_Nb = 1.102
coeff_a5_low_energy_alpha_empirical_formula_Mo = 1.29
coeff_a5_low_energy_alpha_empirical_formula_Tc = 5.792
coeff_a5_low_energy_alpha_empirical_formula_Ru = 5.754
coeff_a5_low_energy_alpha_empirical_formula_Rh = 8.689
coeff_a5_low_energy_alpha_empirical_formula_Pd = 1.289
coeff_a5_low_energy_alpha_empirical_formula_Ag = 2.844
coeff_a5_low_energy_alpha_empirical_formula_Cd = 3.119
coeff_a5_low_energy_alpha_empirical_formula_In = 6.026
coeff_a5_low_energy_alpha_empirical_formula_Sn = [7.767,7.767]
coeff_a5_low_energy_alpha_empirical_formula_Sb = 2.723
coeff_a5_low_energy_alpha_empirical_formula_Te = 3.411
coeff_a5_low_energy_alpha_empirical_formula_I  = 1.51
coeff_a5_low_energy_alpha_empirical_formula_Xe = 0.9126
coeff_a5_low_energy_alpha_empirical_formula_Cs = 2.597
coeff_a5_low_energy_alpha_empirical_formula_Ba = 3.726
coeff_a5_low_energy_alpha_empirical_formula_La = 0.9121
coeff_a5_low_energy_alpha_empirical_formula_Ce = 2.529
coeff_a5_low_energy_alpha_empirical_formula_Pr = 2.444
coeff_a5_low_energy_alpha_empirical_formula_Nd = 2.328
coeff_a5_low_energy_alpha_empirical_formula_Pm = 2.258
coeff_a5_low_energy_alpha_empirical_formula_Sm = 1.997
coeff_a5_low_energy_alpha_empirical_formula_Eu = 1.016
coeff_a5_low_energy_alpha_empirical_formula_Gd = 0.9244
coeff_a5_low_energy_alpha_empirical_formula_Tb = 1.258
coeff_a5_low_energy_alpha_empirical_formula_Dy = 1.055
coeff_a5_low_energy_alpha_empirical_formula_Ho = 1.221
coeff_a5_low_energy_alpha_empirical_formula_Er = 0.849
coeff_a5_low_energy_alpha_empirical_formula_Tm = 1.604
coeff_a5_low_energy_alpha_empirical_formula_Yb = 1.459
coeff_a5_low_energy_alpha_empirical_formula_Lu = 1.66
coeff_a5_low_energy_alpha_empirical_formula_Hf = 1.589
coeff_a5_low_energy_alpha_empirical_formula_Ta = 0.8619
coeff_a5_low_energy_alpha_empirical_formula_W  = 0.8228
coeff_a5_low_energy_alpha_empirical_formula_Re = 1.263
coeff_a5_low_energy_alpha_empirical_formula_Os = 1.191
coeff_a5_low_energy_alpha_empirical_formula_Ir = 1.305
coeff_a5_low_energy_alpha_empirical_formula_Pt = 0.8962
coeff_a5_low_energy_alpha_empirical_formula_Au = 1.05
coeff_a5_low_energy_alpha_empirical_formula_Hg = 1.416
coeff_a5_low_energy_alpha_empirical_formula_Tl = 1.386
coeff_a5_low_energy_alpha_empirical_formula_Pb = 3.224

# from Berger et al. - Stopping powers and ranges for protons and alpha particles - ICRU report 49 (1993) - Table 3.7
# Rk : I imposed : T1 = 30 MeV and T2 = 40 MeV for Lithium    [Li]3,
#                  T1 =  2 MeV and T2 =  3 MeV for Berylium   [Be]4,
#                  T1 =  5 MeV and T2 = 10 MeV for Oxygen     [O]8,
#                  T1 =  5 MeV and T2 = 10 MeV for Fluorine   [F]9,
#                  T1 =  5 MeV and T2 = 10 MeV for Neon       [Ne]10,
#                  T1 = 30 MeV and T2 = 40 MeV for Sodium     [Na]11,
#                  T1 = 10 MeV and T2 = 20 MeV for Zirconium  [Zr]40,
#                  T1 = 10 MeV and T2 = 20 MeV for Niobium    [Nb]41,
#                  T1 = 10 MeV and T2 = 20 MeV for Molybdenum [Mo]42,
#                  T1 = 10 MeV and T2 = 20 MeV for Technetium [Tc]43,
#                  T1 = 10 MeV and T2 = 20 MeV for Ruthenium  [Ru]44,
#                  T1 = 10 MeV and T2 = 20 MeV for Rhodium    [Rh]45,
#                  T1 = 10 MeV and T2 = 20 MeV for Palladium  [Pd]46,
#                  T1 = 10 MeV and T2 = 20 MeV for Silver     [Ag]47,
#                  T1 = 10 MeV and T2 = 20 MeV for Cadmium    [Cd]48,
#                  T1 = 10 MeV and T2 = 20 MeV for Indium     [In]49,
#                  T1 = 10 MeV and T2 = 20 MeV for Tin        [Sn]50, ... until [Pb]82

param_T1_low_energy_alpha_empirical_formula_H  = 1.0
param_T1_low_energy_alpha_empirical_formula_He = 1.5
param_T1_low_energy_alpha_empirical_formula_Li = 30. #1.0
param_T1_low_energy_alpha_empirical_formula_Be = 2.0 #1.0
param_T1_low_energy_alpha_empirical_formula_B  = 1.0
param_T1_low_energy_alpha_empirical_formula_C  = [1.0,1.0,1.0]
param_T1_low_energy_alpha_empirical_formula_N  = 1.0
param_T1_low_energy_alpha_empirical_formula_O  = 5. #1.0
param_T1_low_energy_alpha_empirical_formula_F  = 5.0 #1.0
param_T1_low_energy_alpha_empirical_formula_Ne = 5.0 #1.0
param_T1_low_energy_alpha_empirical_formula_Na = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Mg = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Al = 2.0
param_T1_low_energy_alpha_empirical_formula_Si = 2.0
param_T1_low_energy_alpha_empirical_formula_P  = [30.,30.,30.,30.] #[0.8,0.8,0.8,0.8]
param_T1_low_energy_alpha_empirical_formula_S  = [0.8,0.8,0.8]
param_T1_low_energy_alpha_empirical_formula_Cl = 0.8
param_T1_low_energy_alpha_empirical_formula_Ar = 0.8
param_T1_low_energy_alpha_empirical_formula_K  = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Ca = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Sc = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Ti = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_V  = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Cr = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Mn = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Fe = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Co = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Ni = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Cu = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Zn = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Ga = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_Ge = 30. #2.0
param_T1_low_energy_alpha_empirical_formula_As = 3.0
param_T1_low_energy_alpha_empirical_formula_Se = [3.0,3.0,3.0]
param_T1_low_energy_alpha_empirical_formula_Br = 3.0
param_T1_low_energy_alpha_empirical_formula_Kr = 3.0
param_T1_low_energy_alpha_empirical_formula_Rb = 2.0
param_T1_low_energy_alpha_empirical_formula_Sr = 2.0
param_T1_low_energy_alpha_empirical_formula_Y  = 2.0
param_T1_low_energy_alpha_empirical_formula_Zr = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Nb = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Mo = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Tc = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Ru = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Rh = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Pd = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Ag = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Cd = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_In = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Sn = [10.,10.] #2.0
param_T1_low_energy_alpha_empirical_formula_Sb = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Te = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_I  = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Xe = 10. #2.0
param_T1_low_energy_alpha_empirical_formula_Cs = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Ba = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_La = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Ce = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Pr = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Nd = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Pm = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Sm = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Eu = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Gd = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Tb = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Dy = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Ho = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Er = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Tm = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Yb = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Lu = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Hf = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Ta = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_W  = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Re = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Os = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Ir = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Pt = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Au = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Hg = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Tl = 10. #1.0
param_T1_low_energy_alpha_empirical_formula_Pb = 10. #1.0

param_T2_low_energy_alpha_empirical_formula_H  = 2.0
param_T2_low_energy_alpha_empirical_formula_He = 3.0
param_T2_low_energy_alpha_empirical_formula_Li = 40. #4.0
param_T2_low_energy_alpha_empirical_formula_Be = 4. #4.0
param_T2_low_energy_alpha_empirical_formula_B  = 3.0
param_T2_low_energy_alpha_empirical_formula_C  = [3.0,3.0,3.0]
param_T2_low_energy_alpha_empirical_formula_N  = 3.0
param_T2_low_energy_alpha_empirical_formula_O  = 10. #2.0
param_T2_low_energy_alpha_empirical_formula_F  = 10. #2.0
param_T2_low_energy_alpha_empirical_formula_Ne = 10. #2.0
param_T2_low_energy_alpha_empirical_formula_Na = 50. #4.0
param_T2_low_energy_alpha_empirical_formula_Mg = 50. #4.0
param_T2_low_energy_alpha_empirical_formula_Al = 4.0
param_T2_low_energy_alpha_empirical_formula_Si = 4.0
param_T2_low_energy_alpha_empirical_formula_P  = [50.,50.,50.,50.]#[3.0,3.0,3.0,3.0]
param_T2_low_energy_alpha_empirical_formula_S  = [3.0,3.0,3.0]
param_T2_low_energy_alpha_empirical_formula_Cl = 3.0
param_T2_low_energy_alpha_empirical_formula_Ar = 3.0
param_T2_low_energy_alpha_empirical_formula_K  = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Ca = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Sc = 20. #6.0
param_T2_low_energy_alpha_empirical_formula_Ti = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_V  = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Cr = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Mn = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Fe = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Co = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Ni = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Cu = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Zn = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Ga = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_Ge = 50. #6.0
param_T2_low_energy_alpha_empirical_formula_As = 7.0
param_T2_low_energy_alpha_empirical_formula_Se = [7.0,7.0,7.0]
param_T2_low_energy_alpha_empirical_formula_Br = 7.0
param_T2_low_energy_alpha_empirical_formula_Kr = 7.0
param_T2_low_energy_alpha_empirical_formula_Rb = 5.0
param_T2_low_energy_alpha_empirical_formula_Sr = 5.0
param_T2_low_energy_alpha_empirical_formula_Y  = 5.0
param_T2_low_energy_alpha_empirical_formula_Zr = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Nb = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Mo = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Tc = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Ru = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Rh = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Pd = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Ag = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Cd = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_In = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Sn = [20.,20.] #5.0
param_T2_low_energy_alpha_empirical_formula_Sb = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Te = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_I  = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Xe = 20. #5.0
param_T2_low_energy_alpha_empirical_formula_Cs = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Ba = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_La = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Ce = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Pr = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Nd = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Pm = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Sm = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Eu = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Gd = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Tb = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Dy = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Ho = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Er = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Tm = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Yb = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Lu = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Hf = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Ta = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_W  = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Re = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Os = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Ir = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Pt = 20. #3.0
param_T2_low_energy_alpha_empirical_formula_Au = 20. #2.0
param_T2_low_energy_alpha_empirical_formula_Hg = 20. #2.0
param_T2_low_energy_alpha_empirical_formula_Tl = 20. #2.0
param_T2_low_energy_alpha_empirical_formula_Pb = 20. #2.0
