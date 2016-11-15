import numpy as np
from astropy.table import Table,Column
import os,sys
import time
from matplotlib import pyplot as plt
from scipy.stats import skewnorm
from scipy.interpolate import UnivariateSpline
from astropy.cosmology import FlatLambdaCDM, FlatwCDM
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

setplot = 0 #debugging flag to plot results
bench = 0 #debugging flag to print to screen


class SncosmoSimulation(object):
	def __init__(self,survey_fields=None, simlib_file=None,simlib_obs_sets=None, c_pdf='Normal',\
			 c_params=[-0.05,0.2],x1_pdf='Normal',x1_params=[0.5,1.0],\
			 t0min = 56525.0,t0max=57070.0,NumSN = 500,minNumEpochs = 5,minNumFilters = 3, minSNR = 4.0,\
			cosmo=None,alpha=0.14,beta=3.2,deltaM=0.0,zp_off=[0.0,0.0,0.0,0.0]):
		'''Input: 
			survey_fields:dict, survey fields to generate z dist, values should be [area,zmin,zmax,dndz func,time ]
			simlib_file: str, snana simlib  
			simlib_obs_sets= observation library
			c_pdf and x1_pdf: str, either 'Normal' or 'SkewNormal'
			c_params and x1_params: list, hyperparams for 'Normal' [mean,var] or 'SkewNormal' [mean,var,skew]
			t0min: float
			t0max: float
			NumSN: int, number of SN to simulate and fit
			minNumEpochs: int, minNumFilters: int, minSNR: int,  selection cut: require minNumEpochs in at least minNumFilters with SNR>minSNR 
			cosmo: astropy cosmology
			alpha: float
			beta: float
			deltaM: magnitude offset parameter
			zp_off: list, zp offsets in each band
		'''
		self.t0min=t0min ; self.t0max = t0max
		self.NumSN = NumSN
		self.minNumEpochs = minNumEpochs
		self.minNumFilters = minNumFilters
		self.minSNR = minSNR
		self.alpha = alpha ;self.beta=beta;self.deltaM=deltaM ;self.zp_off = zp_off

		if survey_fields ==None:
			self.survey_fields = self.DES_specific_zdist()
		else:
			self.survey_fields = survey_fields
		if cosmo ==None: 
			self.cosmo=FlatLambdaCDM(H0=70.0, Om0=0.3)
		else:
			self.cosmo = cosmo
		self.totz=self.get_zdist()
		print "Total number of sn to be simulated:",len(self.totz)
		
		start = time.time()
		if simlib_file: 
			self.simlib_meta, self.simlib_obs_sets = self.read_simlib(simlib_file)
		else:
			self.simlib_obs_sets =simlib_obs_sets 
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
		self.model = sncosmo.Model(source='salt2-extended',effects=[dust],effect_names=['mw'],\
		effect_frames=['obs'])
		self.get_parameter_list()
		self.generate_lcs()
		start = time.time()
		self.fit_lcs()
		end = time.time()
		if bench: print "Fitting took", end-start, "secs"
	
			
		if setplot:
			for ii in range(len(self.fit_results)):
				plt.plot(self.simx0[ii],self.fit_results[ii]['parameters'][2],\
				linestyle='None', marker="o")
			plt.savefig('x0.png')
			plt.close()
			for ii in range(len(self.fit_results)):
				plt.plot(self.simx1[ii],self.fit_results[ii]['parameters'][3],\
				linestyle='None', marker="o")
			plt.savefig('x1.png')
			plt.close()
			for ii in range(len(self.fit_results)):
				plt.plot(self.simc[ii],self.fit_results[ii]['parameters'][4],\
				linestyle='None', marker="o")
			plt.savefig('c.png')
			plt.close()
		
		#print self.fit_results[0].keys(), self.fit_results[0]['vparam_names'], self.fit_results[0]['parameters']
		#print self.fit_results[0]['parameters']
		#print self.fit_results[0]['covariance'], self.fit_results[0]['covariance'].shape
		#print self.totz[0],self.simt0[0], self.simx0[0], self.simx1[0], self.simc[0]
		if setplot:
			sncosmo.plot_lc(self.lcs[0][0], model=self.fitted_model[0], errors=self.fit_results[0].errors)
			plt.savefig("lc.png")
		

	def fit_lcs(self):
		self.fit_results=[]
		self.fitted_model=[]
		for jj,lc in enumerate(self.lcs):
    			self.model.set(z=self.totz[jj])
			res, fitted_model = sncosmo.fit_lc(lc[0], self.model,['z','x0', 'x1', 'c','t0'],\
			bounds={'x1':(-3.0, 3.0), 'c':(-0.3,0.3),'z':(0.0,1.2)}, minsnr =self.minSNR)
			self.fit_results.append(res)
			self.fitted_model.append(fitted_model)


	def impose_SNR_cuts(self,lc):
		'''require minNumEpochs in at least minNumFilters with SNR>minSNR'''
		count=0
		bands = ['desg','desr','desi','desz']
		for b in bands:
			if len(lc['time'][(lc['flux']/lc['fluxerr'] > self.minSNR) & (lc['band'] == "desg")])  >=self.minNumEpochs:
				count = count+1
		if count >= self.minNumFilters: 
			return 1
		else:
			return 0
			

	def generate_lcs(self):
		if bench: print "Generating lcs...."
		self.lcs=[]
		start = time.time()
		for p in self.params:
			libid = np.random.choice(self.simlib_obs_sets.keys())
			light_curve = sncosmo.realize_lcs(self.simlib_obs_sets[libid], self.model, [p])
			flag = self.impose_SNR_cuts(light_curve[0])
			if flag:
				for jj, b in enumerate(light_curve[0]['band']):
					if b == "desg": light_curve[0]['zp'][jj] += self.zp_off[0]
					if b == "desr": light_curve[0]['zp'][jj] += self.zp_off[1]
					if b == "desi": light_curve[0]['zp'][jj] += self.zp_off[2]
					if b == "desz": light_curve[0]['zp'][jj] += self.zp_off[3]
				self.lcs.append(light_curve)
		end = time.time()
		if bench:
			print "Simulation took", end-start, "secs"
			print "Number of simulations passing cuts:", len(self.lcs)
		

	def get_parameter_list(self):
		'''Create list of parameter dictionaries '''
		self.params = []
		self.mu=[]
		self.simx0 = np.zeros(len(self.simc))
		for ii,z in enumerate(self.totz):
    			self.model.set(z=z)
			mu = self.cosmo.distmod(z)
			self.mu.append(mu.value)
			mb = mu.value - self.alpha*self.simx1[ii] + self.beta*self.simc[ii] - 19.35 - self.deltaM
			self.model.source.set_peakmag(mb, 'bessellb', 'ab')
    			self.simx0[ii] = self.model.get('x0')
    			p = {'z':z, 't0':self.simt0[ii], 'x0':self.simx0[ii], 'x1':self.simx1[ii] , 'c': self.simc[ii]}
    			self.params.append(p)
	
	def generate_t0(self):
		''' function to generate t0 from uniform dist '''
		if bench: 
			print 'Drawing t0 from uniform dist...'
		self.simt0 = np.random.uniform(self.t0min,self.t0max,size=len(self.totz))

	def generate_random_c_x1(self):
		'''Function to generate random c and x1 from either Normal or skew normal distributions'''
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
		
	@staticmethod
	def read_simlib(libfile):
		'''function to read snana like simlib file'''
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
    		return 2.6E-5*math.pow((1+z),1.5)


	def snana_snrate_hiz(self,z):
    		'''FROM SNANA_MAIN_FILE.input'''
    		return 7.35E-5

	def DES_specific_zdist(self):
		'''create dict specific to DES fields, areas are in deg^2, 
		survey time is specified from GENRANGE_PEAKMJD, zmin, zmax and zcut are specified in SNANA GENRANGE_REDSHIFT & DNDZ 
		'''
        	xarea, carea, earea, sarea = 17.1173, 16.2981,11.6045,12.7980 
        	surveytime = 525 #540 
		zmin, zmax, zcut = 0.05, 1.2, 1.0 
        	DES_fields = {'xlo':[xarea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'xhi':[xarea,zcut,zmax,self.snana_snrate_hiz,surveytime],\
			'clo':[carea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'chi':[carea,zcut,zmax,self.snana_snrate_hiz,surveytime],\
			'elo':[earea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'ehi':[earea,zcut,zmax,self.snana_snrate_hiz,surveytime],\
			'slo':[sarea,zmin,zcut,self.snana_snrate_loz,surveytime],\
			'shi':[sarea,zcut,zmax,self.snana_snrate_hiz,surveytime]}

        	return DES_fields

	def inv_transform_spline(self,z):
		'''Method to create inverse spline to discrete cumulative distribution function 
		to allow drawing random variables.
		Warning: user should check that spline faithfully matches actual cdf.
		'''
		srt_param=np.sort(z)
		cdf = np.array(range(len(z)))/float(len(z))
		#create a spline
		self.spline2_cdf = UnivariateSpline(cdf,srt_param,k=5)


	def draw_fracZ_from_field(self,z,frac):
		num = frac*self.NumSN
		self.inv_transform_spline(z)
		subsampled_z =[]
		for n in range(0,int(num)):
			uni_rvs = np.random.uniform()
			subsampled_z.append(float(self.spline2_cdf(uni_rvs)))
		return subsampled_z
		

	def get_zdist(self):
		sn_z=[]
		for field in self.survey_fields.values():
			z=list(sncosmo.zdist(field[1],field[2],area=field[0],ratefunc=field[3],time=field[4],cosmo = self.cosmo))
			sn_z = np.concatenate((sn_z,z))
		total_z = len(sn_z)
		truncated_z=[]
		for field in self.survey_fields.values():
                        z=list(sncosmo.zdist(field[1],field[2],area=field[0],ratefunc=field[3],time=field[4],cosmo=self.cosmo))
			field_frac = len(z)/float(total_z)
			z_sample = self.draw_fracZ_from_field(z,field_frac)
			truncated_z = np.concatenate((truncated_z,z_sample))

		return truncated_z

if __name__=='__main__':

	#Read simlib file before sampling, obs is list of libIDs which is passed to simulation class each time.
	#This will take a while as we also need to sort all times in all bands so they are monotonically increasing
	libfile = 'DES_DIFFIMG.SIMLIB' #NOTE: THE SIMLIB INDEX STARTS AT 1, NOT 0#
	meta,obs=SncosmoSimulation.read_simlib(libfile)
	for libid in obs.keys():
		 df = obs[libid].to_pandas()
		 df=df.sort_values(by='time')
		 obs[libid] = Table.from_pandas(df)

	#parameters we might want to vary in sampler
	cosmo = FlatwCDM(name='SNLS3+WMAP7', H0=71.58, Om0=0.262, w0=-1.0)
	alpha=0.14
	beta=3.2
	deltaM = 0.0
	zp_offset = [0.0,0.0,0.0,0.0] #griz

	fm_sim = SncosmoSimulation(simlib_obs_sets=obs,cosmo = cosmo,alpha=alpha,beta=beta,deltaM=deltaM,zp_off = zp_offset,NumSN =50)


