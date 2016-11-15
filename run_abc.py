
############## superABC ##################################
# Author: Elise Jennings, elise@fnal.gov
# Code to run ABC sampler with SN data and two metrics: Tripp and Light Curve metric.
##########################################################

import numpy as np
from math import floor
import sys
import subprocess
import os
import re
from time import sleep
from itertools import chain
from write_simfit_files import *
from sncosmo_forwardmodel import *
import sys

from parse_params import *
import astroabc as abc

ROOT_dir = os.path.split(os.path.abspath(__file__))[0]
param_file = os.path.join(ROOT_dir,'params.ini')


__all__ = ['main','data_sim_match','scale_sim','dist_metric_tripp_snana', 'dist_metric_tripp_sncosmo','get_ref_pdf','dist_metric_LC','scale_and_spline']


def data_sim_match(dLC,sLC):
        '''Light Curve metric. Create a dict with keys=SimID_4LCs, value=[[DataID_4LCs,pkmjd-rescale,zrescale]]
        input:
                dLC: instance of lightcurves for the data
                sLC: instance of lightcurves for one sim
        returns:
                sim_dict: dict of key,value = SimID_4LCs, value=[[DataID_4LCs,pkmjd-rescale,zrescale]]
        '''
        sim_dl = sLC.dl_spline

        sim_dict = {}
        for scid in sLC.fitres.keys():

                idx = np.where(np.array(dLC.simlib_libid) == sLC.fitres[scid][2])[0] #index in the cids array
                dcid = int(dLC.simlib_cids[idx])

                # find z rescaling needed due to z diff
                old_z = sLC.fitres[scid][1] ; new_z = dLC.fitres[dcid][1]
                z_rescale = float((sim_dl(old_z)/sim_dl(new_z))**2)

                val = [dLC.fitres[dcid][2],float(dLC.fitres[dcid][0] - sLC.fitres[scid][0]), z_rescale]

                indx_lc = sLC.fitres[scid][3]
                sim_dict[indx_lc]=val
        return sim_dict


def get_ref_pdf():
        '''Light Curve metric. Function which reads pre-defined reference difference PDF
            Input: none
            Returns: average ref diff PDF (normalised)
        '''
	ref = np.loadtxt("ref_pdf.txt")
	return ref
            

def dist_metric_LC(dataLC,simLC):
        ''' Light Curve metric. Function to find the observed difference PDF all bands for each SIM-DATA pair.
            This is compared with the expected distribution from the ref PDF to calculate the chi^2
        Input: dataLC - instance of data LC class
                simLC - instance of sim LC class    
        Returns: TS/real_ddof - chi^2 per degree of freedom
        '''

        print "-->Calculating chi2..."

        if simLC == None:
            return np.infty

        sim_dict = data_sim_match(dataLC,simLC)

        xx = []
        sample_num = 0
        Nbins =400 
        for nn,scid in enumerate(sim_dict.keys()):
            vals  = sim_dict[scid]
            d_mjd,d_f,d_ferr = dataLC.get_obs_lc(vals[0],1)

            s_f = scale_and_spline(scid,vals,d_mjd,simLC)
            for i in range(4):
                xx.append((s_f[i]-d_f[i])/d_ferr[i])
                sample_num = sample_num + len(d_f[i])
        xx = np.fromiter(chain.from_iterable(xx),dtype='float')

        obs_hist,obs_bins = np.histogram(xx, bins=Nbins, range=(-200,200))
        binwidth  = 400./float(Nbins)

        #save ref dist to a file
        #hist,bins = np.histogram(xx, bins=Nbins, range=(-200,200),density=True)
        #data = np.vstack((bins[:-1],hist)).T
        #np.savetxt("reference_distribution.txt",data)

        #load ref dists
        #ref_hist = get_ref_pdf()


        nonzero = np.where(ref_hist > 0.0)[0]
        real_ddof = float(len(nonzero))
        TS = np.sum((obs_hist[nonzero] - ref_hist[nonzero]*sample_num*binwidth)**2/(ref_hist[nonzero]*sample_num*binwidth))

        print "-->...Done"
        return TS/real_ddof  #chi2/dof



def scale_and_spline(scid,val,d_mjd,sLC):
        '''Light Curve metric. Function to obtain a sim SN LC and adjust it to match data pkmjd and redshift in each band.
        Input: scid - id of the sim SN
                val - list of data id, pkmjd correction, z correction
                d_mjd - epochs of data LC
                sLC - instance of sim SN LC class
        Returns: list of lists of splined sim SN fluxes at same epochs as data in each band
        '''
        dcid, mjd_corr,z_corr = val
        s_mjd,s_f_old,s_ferr = sLC.get_obs_lc(int(scid),0)
        splined_fluxes = []
        for i,dx in enumerate(d_mjd):
            s_mjd_new = s_mjd[i] + mjd_corr
            s_f_new = np.array(s_f_old[i])*z_corr
            sim_sp = sLC.create_spline(s_mjd_new,s_f_new,s_ferr[i])
            y = sim_sp(dx)
            splined_fluxes.append(y)
        return splined_fluxes

def dist_metric_tripp_snana(dataLC,simLC):
    """Tripp metric using snana simulations
    """
    if simLC == None:
        return np.infty

    c_d=np.zeros(len(dataLC.fitres.keys()));x1_d=np.zeros(len(dataLC.fitres.keys()));mu=np.zeros(len(dataLC.fitres.keys()))
    mb_d=np.zeros(len(dataLC.fitres.keys()));z_d=np.zeros(len(dataLC.fitres.keys())); mb_err=np.zeros(len(dataLC.fitres.keys()))
    c_err=np.zeros(len(dataLC.fitres.keys())); x1_err=np.zeros(len(dataLC.fitres.keys()))

    for nn,dcid in enumerate(dataLC.fitres.keys()):
            z_d[nn] = dataLC.fitres[dcid][1]
            x1_d[nn] = dataLC.fitres[dcid][3]
            x1_err[nn] = dataLC.fitres[dcid][4]
            c_d[nn] = dataLC.fitres[dcid][5]
            c_err[nn] = dataLC.fitres[dcid][6]
            mb_d[nn] = dataLC.fitres[dcid][7]
            mb_err[nn] = dataLC.fitres[dcid][8]
            mu[nn] = simLC.mu_func(dataLC.fitres[dcid][1])
    err= np.sqrt((mb_err)**2 + (simLC.alpha_param*x1_err)**2+(simLC.beta_param*c_err)**2 + 0.11**2)
    d_1  = 1./np.float(len(mu))*np.sum((mu - (mb_d + simLC.alpha_param*x1_d - simLC.beta_param*c_d + 19.35 +simLC.param_magoff))**2/err**2)

    c_s=np.zeros(len(simLC.fitres.keys()));x1_s=np.zeros(len(simLC.fitres.keys()));mu_s=np.zeros(len(simLC.fitres.keys()))
    mb_s=np.zeros(len(simLC.fitres.keys())); mb_s_err=np.zeros(len(simLC.fitres.keys()))
    c_s_err=np.zeros(len(simLC.fitres.keys())); x1_s_err=np.zeros(len(simLC.fitres.keys()))

    for nn,scid in enumerate(simLC.fitres.keys()):
            x1_s[nn] = simLC.fitres[scid][10]
            x1_s_err[nn] = simLC.fitres[scid][11]
            c_s[nn] = simLC.fitres[scid][8]
            c_s_err[nn] = simLC.fitres[scid][9]
            mb_s[nn] = simLC.fitres[scid][12]
            mb_s_err[nn] = simLC.fitres[scid][13]
            mu_s[nn] = simLC.fitres[scid][6]
    err_s= np.sqrt((mb_s_err)**2 + (simLC.alpha_param*x1_s_err)**2+(simLC.beta_param*c_s_err)**2 + 0.11**2)
    d_2  = 1./np.float(len(mu_s))*np.sum((mu_s - (mb_s + simLC.alpha_param*x1_s - simLC.beta_param*c_s + 19.35+simLC.param_magoff))**2/err_s**2)


    metric = np.abs(d_1-d_2)
    node_id = abc.mp.current_process().pid

    print "Done...", metric, node_id
    return metric

def dist_metric_tripp_sncosmo(data_c,sim_c):
    """Tripp metric using sncosmo simulations
    """
    print "finding metric..."
    if sim_c == None:
        return np.infty

    len_data = len(data_c.fit_results)
    len_sim = len(sim_c.fit_results)

    err=np.zeros(len_data); d_1 = np.zeros(len_data)
    err_s=np.zeros(len_sim); d_2 = np.zeros(len_sim)
	

    for i,lc in enumerate(data_c.fit_results):
	err[i] = (-2.5*np.log10(lc['errors']['x0']))**2 + (sim_c.alpha*lc['errors']['x1'] )**2 \
	+ (sim_c.beta*lc['errors']['c'] )**2 + 0.11**2   
	mu = sim_c.cosmo.distmod(data_c.totz[i])
	d_1[i] = mu.value - \
	( -2.5*np.log10(lc['parameters'][1]) + sim_c.alpha*lc['parameters'][2] - \
	sim_c.beta*lc['parameters'][3] + 19.35 + sim_c.deltaM )

    t1= 1./np.float(len_data)*np.sum((d_1)**2/err)

    for i,lc in enumerate(sim_c.fit_results):
	err_s[i] =  (-2.5*np.log10(lc['errors']['x0']))**2 + (sim_c.alpha*lc['errors']['x1'] )**2 \
	 + (sim_c.beta*lc['errors']['c'] )**2 + 0.11**2  
	d_2[i] = sim_c.mu[i] - \
	( -2.5*np.log10(lc['parameters'][1]) + sim_c.alpha*lc['parameters'][2] - \
	sim_c.beta*lc['parameters'][3] + 19.35 + sim_c.deltaM )

    t2= 1./np.float(len_sim)*np.sum((d_2)**2/err_s)
    print t1,t2

    metric = np.abs(t1-t2)
    node_id = abc.mp.current_process().pid
    print " \t \t Done...", metric, node_id
    return metric

def run_simulation(theta):
        '''This is the main function the abc sampler calls for every particle.
        Given a parameter set theta, a sim is run with these parameters and fit
.
        Input: theta = array of parameters to be varied specified in params.ini
        Returns: instance of LightCurves for a simulation using theta, which is passed to dist_metric
        '''
        print "theta", theta

        #option to impose hard boundary to prior range
        #if np.any([theta[i] <hyperp[i][0] or theta[i] > hyperp[i][1] for i in range(len(theta))]):
        #    return None

        node_id = abc.mp.current_process().pid
        theta_om , theta_w0, theta_alpha, theta_beta, theta_magoff,\
        theta_ZPg, theta_ZPr, theta_ZPi, theta_ZPz  = check_order_params(param_names,theta)
        params = (theta_om,1-theta_om,.7,theta_w0,0.)

	#example to impose physical prior
	if theta_om<0.0 or theta_om >1.0:
		return None

        if nparam == 1 and param_names != 'w0':
                wmode=-1
        else:wmode=0

        if fm_sim == "snana":
            # launching SNANA sim, this requires us to write out a sim and fit file with the parameters needed.
            sim_root_file = sim_dirname+'OUT_ABC_'+str(node_id)+'/EJ_ABC_'+str(node_id)+'/FITOPT000.ROOT'
            summary_file = sim_dirname+"/SIMLOGS_"+str(node_id)+"/SIMJOB_SUMMARY.LOG"
            merge_file = sim_dirname+"/OUT_ABC_"+str(node_id)+"/MERGE.LOG"

            write_simfile(theta_om,1.-theta_om,theta_w0,theta_beta,theta_alpha, theta_magoff,node_id)
        
            sim_ok = launch_sim(node_id,summary_file)
            if not sim_ok:
                return None
            write_nml(node_id, theta_ZPg, theta_ZPr, theta_ZPi, theta_ZPz)
            launch_fit(node_id,sim_root_file)
            fit_done = False
            while not fit_done:
                fit_done = check_fit(node_id,merge_file)
            sleep(30.) #sometimes we have to wait for the ROOT file to be written out fully...
            simLC = LightCurves("sim",fromfiles=True,ROOTfile=sim_root_file,cosmo_params = params,wmode=wmode,beta=theta_beta,alph=theta_alpha)
            simLC.param_magoff = theta_magoff
            return simLC
        
        elif fm_sim == "sncosmo":

		sim_cosmo = FlatwCDM(name='SNLS3+WMAP7', H0=70., Om0=theta_om, w0=theta_w0)
		sim_alpha= theta_alpha
		sim_beta= theta_beta
		sim_deltaM= theta_magoff
		sim_zp_offset = [theta_ZPg,theta_ZPr,theta_ZPi,theta_ZPz] #griz
		sim_c = SncosmoSimulation(simlib_obs_sets=obs,cosmo = sim_cosmo,alpha=sim_alpha,beta=sim_beta,\
		deltaM=sim_deltaM,zp_off = sim_zp_offset,NumSN =200)

            	return sim_c

        else:
            print "Please specify the forward model simulation..."
            return 0
            

	

def main():
        """
        Function to initialize instance of the sampler, specify the forward model simulation and begin sampling.
        """
	if fm_sim == "snana":
        	sampler = abc.ABC_class(nparam,npart,dataLC,tlevels,niter,prior,**prop)
	elif fm_sim == "sncosmo":
        	sampler = abc.ABC_class(nparam,npart,data,tlevels,niter,prior,**prop)
        model_sim = run_simulation

        sampler.sample(model_sim)

	
if __name__ == "__main__":

        pconfig = ParseParams(param_file)
        if pconfig.verbose: print "parameters in file: ",  vars(pconfig)        

        sim_dirname = pconfig.sim_dirname
        ref_pdf_mode = pconfig.ref_pdf_mode
        ref_pdf_fname ="None" 
        Rfile = pconfig.rootfile 
        simlib = pconfig.simlib
        fm_sim = pconfig.fm_sim

        print "\n"
        print "Reading data files..."

	if fm_sim == "snana":
		from light_curves import *
        	dataparams=(0.30,0.70,.7,-1,0.)  #Omega_m,Omega_l,h0,w0 #This isn't used for data...dummy input to LightCurves class
        	dataLC = LightCurves("data",fromfiles=True,ROOTfile=Rfile,cosmo_params = dataparams,simlib=simlib)
	elif fm_sim == "sncosmo":
		#read from file here or generate mock dataset (see sncosmo documentation for data structure
		#example here generates mock data
		meta,obs=SncosmoSimulation.read_simlib(simlib)
		for libid in obs.keys():
			df = obs[libid].to_pandas()
			df=df.sort_values(by='time')
			obs[libid] = Table.from_pandas(df)
		#obs is now global and we don't need to read it again

		#parameters for mock data
		cosmo = FlatwCDM(name='SNLS3+WMAP7', H0=70., Om0=0.3, w0=-1.0)
		alpha=0.14
		beta=3.2
		deltaM = 0.0
		zp_offset = [0.0,0.0,0.0,0.0] #griz

        	data = SncosmoSimulation(simlib_obs_sets=obs,cosmo = cosmo,alpha=alpha,beta=beta,deltaM=deltaM,\
		zp_off = zp_offset,NumSN =200)

        print "...Done"

        if pconfig.dfunc=="dist_metric_LC":ref_hist = get_ref_pdf()


        #########
        #Parameters of abc run, these are global
        ########
        nparam =pconfig.nparam
        npart = pconfig.npart
        niter = pconfig.niter
        tlevels = pconfig.tlevels
        out_prefix = pconfig.outfile
        re_prefix = pconfig.restart
        outF = out_prefix+"_"+str(npart)+"part_"+str(niter)+"iter_"+str(nparam)+"nparam.txt"
        restartF = re_prefix+"_abc_"+str(npart)+"part_"+str(niter)+"iter_"+str(nparam)+"nparam.txt"

        if pconfig.num_proc == 1: 
            num_proc = None
        else:
            num_proc = pconfig.num_proc

	dist_method = globals()[pconfig.dfunc]

        prop={'tol_type':pconfig.tol_type,"verbose":pconfig.verbose,'adapt_t':pconfig.adapt_t,'threshold':pconfig.threshold,
        'pert_kernel':pconfig.pert_kernel, 'dist_type': pconfig.dist_type,'dfunc':dist_method,
        'outfile':outF,'num_proc':num_proc, 'from_restart':pconfig.from_restart,'restart':restartF,'mpi':pconfig.mpi,'mp':pconfig.mp}


        ###########
        param_names = np.array(pconfig.pnames)
        priorname  = pconfig.priors
        hyperp = pconfig.hyperp
        prior = zip(priorname,hyperp)

	main()
