#!/usr/bin/env python
import numpy as np
#import random as r
from astropy.table import Table,Column
import os,sys
#from astropy.cosmology import wCDM
#import math
from matplotlib import pyplot as plt
#import scipy as sp
from scipy.stats import skewnorm
import sncosmo
import math
from astropy.io import ascii

##LC SIMULATION AND FITTING USING SNCOSMO
##DESIGN FOR FORWARD-MODEL SIMULATIONS IN
##superABC & BAMBIS

##CODE APPLICATION WRITTEN BY E. JENNINGS & R. WOLF
##SNCOSMO DEVELOPED BY K. BARBARY

##IMPORTANT NOTES
##1. MAKE SURE YOU HAVE SNCOSMO INSTALLED!

##GENERAL NOTES
##1. R.W. LIKES TO PLOT AS SHE WRITES CODE.
##   IF YOU WANT TO SEE PLOTS, SET SETPLOT = 1.
##2. R.W. HAS PUT BENCHMARK PRINT STATEMENTS
##   TO SEE HOW LONG EACH PART TAKES.
##   IF YOU WANT TO SEE PRINTS, SET BENCH = 1.

#########################

setplot = 1
bench = 1

##FUNCTIONS##

###PART 1: THE REDSHIFT DISTRIBUTION###
def snana_snrate_loz(z):
    #FROM SNANA_MAIN_FILE.input & SNANA MANUAL P.45
    return 2.6E-5*math.pow((1+z),1.5)
def snana_snrate_hiz(z):
    #FROM SNANA_MAIN_FILE.input & SNANA MANUAL P.45
    return 7.35E-5

##PART 1: THE REDSHIFT DISTRIBUTION
##---------------------------------
##FROM THE SIMLIB, WE CAN GET THE AREAS OF THE DES FIELDS
##THIS SHOULD PROBABLY BE USER INPUT, BUT I'M JUST GOING TO
##HARD CODE THIS FOR NOW

print 'DRAWING THE REDSHIFT DISTRIBUTION'

xarea, carea, earea, sarea = 17.1173, 16.2981,11.6045,12.7980 #in deg^2
surveytime = 540 #getting this from SNANA_MAIN_FILE.input GENRANGE_PEAKMJD
zmin, zmax, zcut = 0.05, 1.2, 1.0 #from SNANA_MAIN_FILE.input GENRANGE_REDSHIFT & DNDZ

#WE NEED A LOW AND HIGH Z DISTRIBUTION IN EACH OF THE FOUR FIELD GROUPS
x_loz = sncosmo.zdist(zmin,zcut,area=xarea,ratefunc=snana_snrate_loz,time=surveytime)
x_hiz = sncosmo.zdist(zcut,zmax,area=xarea,ratefunc=snana_snrate_hiz,time=surveytime)
c_loz = sncosmo.zdist(zmin,zcut,area=carea,ratefunc=snana_snrate_loz,time=surveytime)
c_hiz = sncosmo.zdist(zcut,zmax,area=carea,ratefunc=snana_snrate_hiz,time=surveytime)
e_loz = sncosmo.zdist(zmin,zcut,area=earea,ratefunc=snana_snrate_loz,time=surveytime)
e_hiz = sncosmo.zdist(zcut,zmax,area=earea,ratefunc=snana_snrate_hiz,time=surveytime)
s_loz = sncosmo.zdist(zmin,zcut,area=sarea,ratefunc=snana_snrate_loz,time=surveytime)
s_hiz = sncosmo.zdist(zcut,zmax,area=sarea,ratefunc=snana_snrate_hiz,time=surveytime)

totz = np.concatenate((list(x_loz),list(x_hiz),list(c_loz),list(c_hiz),list(e_loz),list(e_hiz),list(s_loz),list(s_hiz)))

#testz1 = list(sncosmo.zdist(0,1.0,ratefunc=snana_snrate_loz))
#testz2 = list(sncosmo.zdist(1.0,2.0,ratefunc=snana_snrate_hiz))
#testz = np.concatenate((testz1,testz2))

if(setplot == 1):

    #COMPARING THE SNANA TEST DISTRIBUTION TO THE SNCOSMO DISTRIBUTION
    snana_sim_file = 'SNANA_ALL.DUMP'
    snana = ascii.read(snana_sim_file)
    snana_z = snana['GENZ']
    zbins = np.arange(0,1.3,0.1)
    print 'Size of SNANA ALL:', len(snana_z)
    print 'Size of SNCOSMO ALL:', len(totz)
    plt.hist(totz,normed=True,label='SNCOSMO',bins=zbins)
    plt.hist(snana_z,normed=True,label='SNANA',color='red',histtype='step',lw=2,bins=zbins)
    plt.xlabel('Redshift')
    plt.ylabel('Normed PDF')
    plt.legend(loc='upper left',frameon=False)
    plt.title('Redshift Distribution Before Cuts')
    #plt.hist(testz,normed=True)
    plt.show()

##PART 2: READING IN THE SIMLIB
##---------------------------------
print 'READING THE SNANA SIMLIB. THIS MIGHT TAKE A FEW MINUTES - BE PATIENT!'

snana_simlib_file = 'DES_DIFFIMG.SIMLIB'

snana_simlib_meta, snana_simlib_obs_sets = sncosmo.read_snana_simlib(snana_simlib_file)
#NOTE: THE SIMLIB INDEX STARTS AT 1, NOT 0#

#print snana_simlib_obs_sets.keys()
#print snana_simlib_obs_sets[1].meta
#print snana_simlib_obs_sets[1].colnames

##PART 3: GENERATING THE COLOR
## AND STRETCH DISTRIBUTIONS
##---------------------------------
print 'DRAWING THE COLOR AND STRETCH'

## LET'S ASSUME WE JUST WANT TO DRAW COLOR AND STRETCH
## FROM INDEPENDENT NORMAL DISTRIBUTIONS
simxmean,simxscale = 0.5,1.2
simcmean,simcscale = -0.05, 0.15

simx1 = np.random.normal(loc=simxmean,scale=simxscale,size=len(totz))
simc = np.random.normal(loc=simcmean,scale=simcscale,size=len(totz))

##PART 4: DRAW t0 FROM UNIFORM DISTRIBUTION
##-----------------------------------------
print 'DRAWING t0'
simt0 = np.random.uniform(56525.0,57070,size=len(totz)) #SHOULD AUTOMATE THIS

##PART 5: REALIZE LIGHT CURVES
##----------------------------
print 'REALIZING LIGHT CURVES'
params =[]
for i in xrange(len(totz)):
    p  = {'z':totz[i],'x0':1.0E-5+np.random.normal(loc=0,scale=1E-5), 'x1':simx1[i],'c':simc[i],'t0':simt0[i]}
    #*NOTE: THIS IS A HACK FOR x0 - WE SHOULD ACTUALLY CHAT ABOUT THIS#
    params.append(p)
#for p in params:
#    print p
model = sncosmo.Model(source='salt2')


