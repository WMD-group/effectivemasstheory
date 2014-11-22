#! /usr/bin/env python

"""Calculate some simple semiconductor properties from effective mass theory"""

################################################################################
# Aron Walsh 2014                                                              #
################################################################################

import math as m
import scipy.constants as sc
from numpy import linspace
#import matplotlib.pyplot as plt
from optparse import OptionParser

######################## Set up optional arguments #############################
parser = OptionParser()
parser.add_option("-c", "--electron-effective-mass",
                  action="store", type="float", dest="e", default=0.15,
                  help="Average electron (conduction band) effective mass")
parser.add_option("-v", "--hole-effective-mass",
                  action="store", type="float", dest="h", default=0.14,
                  help="Average hole (valence band) effective mass") 
parser.add_option("-s", "--static-dielectric",
                  action="store", type="float", dest="d0", default=72.01,
                  help="Static (low-frequency) dielectric constant")     
parser.add_option("-o", "--optical-dielectric",
                  action="store", type="float", dest="d1", default=4.96,
                  help="Optical (high-frequency) dielectric constant")    
#########################defaults for CsPbBr3 (PBE0-37.5)#########################        
#parser.add_option("-p", "--optical-phonon",
#                  action="store", type="float", dest="lo", default=9.3, 
#                  help="Optical (polaron active) phonon in THz")
                  
### Further options go here ###
(options,args) = parser.parse_args()

########################### Begin main program #################################
print "*A program to calculate simple semiconductor properties from effective mass theory* \n"
# See, e.g. Fundamentals of Semiconductors, Yu and Cardona

print "Aron Walsh (University of Bath) \nDate last edited: 2/10/2014 \n"

# Get electron effective mass
if options.e ==0:
    e = raw_input("What is the electron effective mass (e.g. 0.3 me)?")
    e = float(e)
else:
    e = options.e

# Get hole effective mass
if options.h ==0:
    h = raw_input("What is the hole effective mass (e.g. 0.3 me)?")
    h = float(h)
else:
    h = options.h
    
# Get static (low frequency) dielectric constant
if options.d0 ==0:
    d0 = raw_input("What is the static dielectric constant (e.g. 10)?")
    d0 = float(d0)
else:
    d0 = options.d0

# Get optical (high frequency) dielectric constant
if options.d1 ==0:
    d1 = raw_input("What is the optical dielectric constant (e.g. 5)?")
    d1 = float(d1)
else:
    d1 = options.d1

# Get optical phonon frequency
#if options.lo ==0:
#    lo = raw_input("What is the optical phonon frequency (e.g. 1 THz)?")
#    lo = float(lo)
#else:
#    lo = options.lo

#
# Calculate properties
#

# Reduced effective mass
    mass=((e*h)/(e+h))
    diel=(1/d1-1/d0)
    print ("*Effective mass \nHole mass: " + str(h) + " me")
    print ("Electron mass: " + str(e) + " me")
    print ("Reduced mass: " + str(mass) + " me \n")

# Exciton Bohr radius
    radius_bohr=(d0/mass)
    radius_bohr_h=(d0/h)
    radius_bohr_e=(d0/e)
    radius=(d0/mass)*0.529177249
    radius_h=(d0/h)*0.529177249
    radius_e=(d0/e)*0.529177249
    print ("*Shallow defects \nAcceptor radius: " + str(radius_h) + " A")
    print ("Donor radius: " + str(radius_e) + " A \n")
    
# (Static) Exciton binding energy
    binding=((1/(d0*radius_bohr))*(13.605698066*1000))
    print ("*Mott-Wannier Analysis \nThermal exciton radius: " + str(radius) + " A")
    print ("Thermal exciton binding energy: " + str(binding) + " meV")    
     
# (Optical) Exciton binding energy
    radius_bohr_o=(d1/mass)
    radius_o=(d1/mass)*0.529177249
    binding_o_ryd=1/(d1*radius_bohr_o)
    binding_o=binding_o_ryd*13.605698066*1000
    print ("\nOptical exciton radius: " + str(radius_o) + " A")
    print ("Optical exciton binding energy: " + str(binding_o) + " meV")        
     
# Carrier polaron radius
    # From Mott (1968)
    radius_bh=(2/(h*diel))*0.529177249
    print ("\nHole (band) polaron radius: " + str(radius_bh) + " A")
    
    radius_be=(2/(e*diel))*0.529177249
    print ("Electron (band) polaron radius: " + str(radius_be) + " A \n")

# Quantum dot properties
    print ("*Quantum Dot")
    
    confine=radius_o*(sc.pi*sc.pi)/3.6
    print ("Confinement radius " + str(confine/10) + " nm")
    
    radius_qd=2.5 #nm
    radius_qd_bohr=radius_qd*18.8971616463
    #radius_ratio=radius_bohr_o/radius_qd_bohr
    #change in band gap (spherical confinement + coulomb attraction + rydberg correction)
    delta_e_ryd=(sc.pi*sc.pi)/(2*mass*radius_qd_bohr*radius_qd_bohr)#-(1.786/(d1*radius_qd_bohr))#-(0.248*binding_o_ryd)
    delta_e=delta_e_ryd*13.605698066*1000
    print ("r=2.5nm band gap correction " + str(delta_e) + " meV")
    
    radius_qd=6 #nm
    radius_qd_bohr=radius_qd*18.8971616463
    #radius_ratio=radius_bohr_o/radius_qd_bohr
    #change in band gap (spherical confinement + coulomb attraction + rydberg correction)
    delta_e_ryd=(sc.pi*sc.pi)/(2*mass*radius_qd_bohr*radius_qd_bohr)#-(1.786/(d1*radius_qd_bohr))#-(0.248*binding_o_ryd)
    delta_e=delta_e_ryd*13.605698066*1000
    print ("r=6nm band gap correction " + str(delta_e) + " meV\n")
    
# Frohlich (lage polaron) properties
    # Speed of light in atomic units
    #   c=1/sc.alpha
    # LO frequency (from THz -> Ry)
#    freq=lo*0.0003039659692
    # Small polaron coupling constant 
#    h_alpha=diel*m.sqrt(h/(2*freq))
#    e_alpha=diel*m.sqrt(e/(2*freq))
    # Small polaron mass (Feynman)
    #   h_pol=h*(1+h_alpha/6)
#    h_pol=h*((1-0.0008*h_alpha*h_alpha)/(1-h_alpha/6+0.0034*h_alpha*h_alpha))
#    radius_bhp=(2/(h_pol*diel))*0.529177249
    #    e_pol=e*(1+e_alpha/6)
#    e_pol=e*((1-0.0008*e_alpha*e_alpha)/(1-e_alpha/6+0.0034*e_alpha*e_alpha))
#    radius_bep=(2/(e_pol*diel))*0.529177249
#    print ("*Hole Polarons \nFrohlich coupling constant: " + str(h_alpha))
#    print ("Effective polaron mass: " + str(h_pol) + " me")
#    print ("Polaron radius: " + str(radius_bhp) + " A \n")
#    print ("*Electron Polarons \nFrohlich coupling constant: " + str(e_alpha))
#    print ("Effective polaron mass: " + str(e_pol) + " me")
#    print ("Polaron radius: " + str(radius_bep) + " A \n")

# Mott transition 

# Exciton transition ~ 1/exciton volume (Optical properties of Solids - Mark Fox)
#    mott=(((0.26/radius_bohr)**3)*(188971616.463**3))
    mott=((1/(4/3*sc.pi*(radius_bohr**3)))*(188971616.463**3))
    print ("*Mott criterion (carrier concentrations) \nExciton: " + str(mott) + " cm-3")

# Mott transition (holes)
    mott=(((0.26/radius_bohr_h)**3)*(188971616.463**3))
    print ("Holes: " + str(mott) + " cm-3")
    
# Mott transition (electrons)
    mott=(((0.26/radius_bohr_e)**3)*(188971616.463**3))
    print ("Electrons: " + str(mott) + " cm-3")
    
# Note that the value of 0.26 for the Mott Criteron is taken from:
# "Universality aspects of the metal-nonmetal transition in condensed media"
# Edwards and Seinko, PRB 17, 2575 (1978) 