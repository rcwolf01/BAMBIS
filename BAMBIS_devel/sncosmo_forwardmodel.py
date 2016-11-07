import numpy as np
from astropy.table import Table,Column
import os,sys
import time
from matplotlib import pyplot as plt
from scipy.stats import skewnorm
try:
	import sncosmo
except ImportError:
	print "Please install sncosmo to use this class"
	sys.exit(0)
import math
from astropy.io import ascii

# Class for LC simulation and fitting using sncosmo
#Written for use in superABC and BAMBIS

##CODE APPLICATION WRITTEN BY E. JENNINGS & R. WOLF
##SNCOSMO DEVELOPED BY K. BARBARY

setplot = 1 #debugging flag to plot results
bench = 1 #debugging flag to print to screen


class SncosmoSimulation(object):
	def __init__(self,survey_fields=None, simlib_file=None, c_pdf='Normal',\
			 c_params=[0.5,1.2],x1_pdf='Normal',x1_params=[-0.05,0.15],\
			 t0min = 56525.0,t0max=57070.0):
		'''Input: 
			survey_fields:dict, survey fields to generate z dist, values should be [area,zmin,zmax,dndz func,time ]
			simlib_file: str, snana simlib  
			c_pdf and x1_pdf: str, either 'Normal' or 'SkewNormal'
			c_params and x1_params: list, hyperparams for 'Normal' [mean,var] or 'SkewNormal' [mean,var,skew]
			t0min: float
			t0max: float
		'''
		self.t0min=t0min ; self.t0max = t0max

		if survey_fields ==None:
			self.survey_fields = self.DES_specific_zdist()
		else:
			self.survey_fields = survey_fields

		self.totz=self.get_zdist()

		start = time.time()
		self.simlib_meta, self.simlib_obs_sets = self.read_simlib(simlib_file)
		end = time.time()
		if bench: print "Read simlib file in ", end-start, "secs"
	
	
		self.c_pdf =c_pdf ; self.c_params=c_params
		self.x1_pdf =x1_pdf ; self.x1_params=x1_params

		start = time.time()
		self.generate_random_c_x1()
		end = time.time()
		if bench: print "c-x1 generated in ", end-start, "secs"
		self.generate_t0()

		#EJ: Note the model needs to be downloaded at this point - is there 
		#anyway we can avoid this for users who are not online?
		dust = sncosmo.CCM89Dust()
		self.model = sncosmo.Model(source='salt2-extended',effects=[dust, dust],effect_names=['host', 'mw'],\
		effect_frames=['rest', 'obs'])
		self.get_parameter_list()
		self.generate_lcs()
		start = time.time()
		self.fit_lcs()
		end = time.time()
		if bench: print "Fitting took", end-start, "secs"

	def fit_lcs(self):
		self.fit_results=[]
		for lc in self.lcs:
			res, fitted_model = sncosmo.fit_lc(lc[0], self.model,['x0', 'x1', 'c'],\
			bounds={'x1':(-3.0, 3.0), 'c':(-0.3,0.3)}, minsnr =3.0)
			self.fit_results.append(res)
		print res.keys()
		print res['param_names']
		print res['parameters']
		print res['covariance']
		print len(self.fit_results)
		

	def generate_lcs(self):
		if bench: print "Generating lcs...."
		self.lcs=[]
		start = time.time()
		for p in self.params:
			libid = 2498 #EJ: Nov 7th right now just generate all with same libid. need to sort times so monotonic so we can use any libid
			#libid = np.random.choice(self.simlib_obs_sets.keys())
			self.lcs.append(sncosmo.realize_lcs(self.simlib_obs_sets[libid], self.model, [p]))
		end = time.time()
		print "Time", end-start
		print type(self.lcs[0][0])
		

	def get_parameter_list(self):
		'''Create list of parameter dictionaries '''
		self.params = []
		for ii,z in enumerate(self.totz):
    			mabs = np.random.normal(-19.3, 0.3)
    			self.model.set(z=z)
    			self.model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    			self.simx0 = self.model.get('x0')
    			p = {'z':z, 't0':self.simt0[ii], 'x0':self.simx0, 'x1':self.simx1[ii] , 'c': self.simc[ii]}
    			self.params.append(p)
	
	def generate_t0(self):
		''' function to generate t0 from uniform dist '''
		if bench: 
			print 'Drawing t0 from uniform dist...'
		self.simt0 = np.random.uniform(self.t0min,self.t0max,size=len(self.totz))

	def generate_random_c_x1(self):
		'''Function to generate random c and x1 from either Normal or skew normal distributions
		'''
		if bench:
			print 'Drawing color and stretch rvs...'
		if self.c_pdf == 'Normal':
			self.simc = np.random.normal(loc=self.c_params[0],scale=self.c_params[1],size=len(self.totz))
		elif self.c_pdf == 'SkewNormal':
			self.simc = skewnorm.rvs(self.c_params[2],loc=self.c_params[0],scale=self.c_params[1],size=len(self.totz))
		if self.x1_pdf == 'Normal':
			self.simx1 = np.random.normal(loc=self.x1_params[0],scale=self.x1_params[1],size=len(self.totz))
		elif self.x1_pdf == 'SkewNormal':
			self.simx1 = skewnorm.rvs(self.x1_params[2],loc=self.x1_params[0],scale=self.x1_params[1],size=len(self.totz))
		

	def read_simlib(self,libfile):
		'''function to read snana like simlib file
		'''
		if libfile==None:
			libfile = 'DES_DIFFIMG.SIMLIB' #NOTE: THE SIMLIB INDEX STARTS AT 1, NOT 0#
		if bench:
			print "Reading the simlib file. This might take a few minutes..."	
		meta,obs= sncosmo.read_snana_simlib(libfile)
		for id in obs.keys():
			#obs[id].rename_column('FLT', 'band')
			obs[id]['band'] = ['desg']*len(obs[id]['FLT'])
			for ii,val in enumerate(obs[id]['FLT']):
				if val =='g': obs[id]['band'][ii] = 'desg'
				if val =='r': obs[id]['band'][ii] = 'desr'
				if val =='i': obs[id]['band'][ii] = 'desi'
				if val =='z': obs[id]['band'][ii] = 'desz'
			obs[id].rename_column('MJD', 'time')
			obs[id].rename_column('CCD_GAIN', 'gain') 
			obs[id].rename_column('SKYSIG', 'skynoise') 
			obs[id].rename_column('ZPTAVG', 'zp') 
			obs[id]['zpsys'] = ['ab']*len(obs[id]['gain'])
		return meta,obs

	def snana_snrate_loz(self,z):
    		'''FROM SNANA_MAIN_FILE.input'''
    		#return 2.6E-5*math.pow((1+z),1.5)
    		return 0.8E-3*math.pow((1+z),0.0)


	def snana_snrate_hiz(self,z):
    		'''FROM SNANA_MAIN_FILE.input'''
    		#return 7.35E-5
    		return 2.6E-5*math.pow((1+z),1.5)

	def DES_specific_zdist(self):
		'''create dict specific to DES fields, areas are in deg^2, 
		survey time is specified from GENRANGE_PEAKMJD, zmin, zmax and zcut are specified in SNANA GENRANGE_REDSHIFT & DNDZ 
		'''
        	xarea, carea, earea, sarea = 17.1173, 16.2981,11.6045,12.7980 
        	surveytime = 525 #540 
		zmin, zmax, zcut = 0.0, 1., 0.08 
        	DES_fields = {'xlo':[xarea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'xhi':[xarea,zcut,zmax,self.snana_snrate_hiz,surveytime],\
			'clo':[carea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'chi':[carea,zcut,zmax,self.snana_snrate_hiz,surveytime],\
			'elo':[earea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'ehi':[earea,zcut,zmax,self.snana_snrate_hiz,surveytime],\
			'slo':[sarea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'shi':[sarea,zcut,zmax,self.snana_snrate_hiz,surveytime]}

        	return DES_fields

	def get_zdist(self):
		sn_z=[]
		for field in self.survey_fields.values():
			z=list(sncosmo.zdist(field[1],field[2],area=field[0],ratefunc=field[3],time=field[4]))
			sn_z = np.concatenate((sn_z,z))

		if setplot:
		#COMPARING THE SNANA TEST DISTRIBUTION TO THE SNCOSMO DISTRIBUTION
    			snana_sim_file = 'SNANA_ALL.DUMP'
    			snana = ascii.read(snana_sim_file)
    			snana_z = snana['GENZ']
    			zbins = np.arange(0,1.3,0.1)
    			print 'Size of SNANA ALL:', len(snana_z)
    			print 'Size of SNCOSMO ALL:', len(sn_z)
    			plt.hist(sn_z,normed=True,label='SNCOSMO',bins=zbins)
    			plt.hist(snana_z,normed=True,label='SNANA',color='red',histtype='step',lw=2,bins=zbins)
    			plt.xlabel('Redshift')
    			plt.ylabel('Normed PDF')
    			plt.legend(loc='upper left',frameon=False)
    			plt.title('Redshift Distribution Before Cuts')
    			plt.savefig('z_dist.png')

		return sn_z

if __name__=='__main__':
	fm_sim = SncosmoSimulation()

