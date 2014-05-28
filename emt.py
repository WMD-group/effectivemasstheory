#! /usr/bin/env python

"""Calculate some simple semiconductor properties from effective mass theory"""

################################################################################
# Aron Walsh 2014                                                              #
################################################################################

from numpy import linspace
#import matplotlib.pyplot as plt
from optparse import OptionParser

######################## Set up optional arguments #############################
parser = OptionParser()
parser.add_option("-c", "--electron-effective-mass",
                  action="store", type="float", dest="e", default=0.15,
                  help="Average electron (conduction band) effective mass")
parser.add_option("-v", "--hole-effective-mass",
                  action="store", type="float", dest="h", default=0.12,
                  help="Average hole (valence band) effective mass") 
parser.add_option("-s", "--static-dielectric",
                  action="store", type="float", dest="d0", default=24.1,
                  help="Static (low-frequency) dielectric constant")     
parser.add_option("-o", "--optical-dielectric",
                  action="store", type="float", dest="d1", default=4.5,
                  help="Optical (high-frequency) dielectric constant")           
### Further options go here ###
(options,args) = parser.parse_args()

########################### Begin main program #################################
print "A program to calculate simple semiconductor properties from effective mass theory"
# See, e.g. Fundamentals of Semiconductors, Yu and Cardona

print "Aron Walsh (University of Bath) \nDate last edited: 25/05/2014 \n"

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

#
# Calculate properties
#

# Reduced effective mass
    mass=((e*h)/(e+h))
    diel=(1/d1-1/d0)
    print ("Reduced effective mass: " + str(mass) + " me \n")

# Exciton Bohr radius
    radius_bohr=(d0/mass)
    radius_bohr_h=(d0/h)
    radius_bohr_e=(d0/e)
    radius=(d0/mass)*0.529177249
    radius_h=(d0/h)*0.529177249
    radius_e=(d0/e)*0.529177249
    print ("Mott-Wannier exciton radius: " + str(radius) + " A")
    print ("Shallow acceptor radius: " + str(radius_h) + " A")
    print ("Shallow donor radius: " + str(radius_e) + " A \n")
    
# Carrier polaron radius
# TO CHECK: Formula has 8*pi --> 2* Prefactor Used
    radius_bh=(2/(h*diel))*0.529177249
    print ("Hole polaron radius: " + str(radius_bh) + " A")
    
    radius_be=(2/(e*diel))*0.529177249
    print ("Electron polaron radius: " + str(radius_be) + " A \n")
    
# Exciton binding energy
    binding=((-1/(2*d0*radius_bohr))*(13.605698066*1000))
    print ("Exciton binding energy: " + str(binding) + " meV")
 
# Note that the value of 0.26 for the Mott Criteron is taken from:
# "Universality aspects of the metal-nonmetal transition in condensed media"
# Edwards and Seinko, PRB 17, 2575 (1978) 
    
# Mott transition (both carriers)
    mott=(((0.26/radius_bohr)**3)*(188971616.463**3))
    print ("Mott criterion (reduced): " + str(mott) + " cm-3")

# Mott transition (holes)
    mott=(((0.26/radius_bohr_h)**3)*(188971616.463**3))
    print ("Mott criterion (holes): " + str(mott) + " cm-3")
    
# Mott transition (electrons)
    mott=(((0.26/radius_bohr_e)**3)*(188971616.463**3))
    print ("Mott criterion (electrons): " + str(mott) + " cm-3")