import numpy as np
from math import floor
import sys
import subprocess
import os
import re
from time import sleep
from itertools import chain
import sys

#EJ: some utility functions for checking sim outputs and writing sim files

def write_simfile(om,ol,w0,beta,alph,magoff,nid):

    '''Function to re-write an SNANA sim file given parameter inputs. New file is labelled with processor ID.
    Input: om=Omega_m, ol=Omega_lambda, w0=w0, beta=color parameter, alph = stretch parameter, nid = processor ID
    '''
    om = round(om,5)
    ol = round(ol,5)
    w0 = round(w0,5)
    beta = round(beta,5)
    alph = round(alph,5)
    magoff = round(magoff,5)
    Filename = 'ONESIM_MASTER.INPUT'
    SIMLIBlines = []
    with open(Filename, 'r') as file:
        SIMLIBlines=file.read().splitlines()
    SIMLIBlines = np.array( SIMLIBlines)

    SIMLIBlines[4] = "GENVERSION: EJ_ABC"+"_"+str(nid)
    newstr ="GENOPT: OMEGA_MATTER "+str(om)+" OMEGA_LAMBDA "+str(ol)+" W0_LAMBDA "+str(w0)+" GENBETA_SALT2 "+str(beta)+" GENALPHA_SALT2 "+str(alph)+" GENMAG_OFF_GLOBAL "+str(magoff)
    SIMLIBlines[21] = "GENPREFIX:   "+str(nid)+"          # prefix of all data filenames"

    Outname = Filename+'_'+str(nid)
    with open(Outname, 'w') as file:
        for l in range(len(SIMLIBlines)):
            if l ==5:
                file.write("%s\n"% newstr)
            else:
                file.write("%s\n" % SIMLIBlines[l])


def write_nml(nid,ZPoff0,ZPoff1,ZPoff2,ZPoff3):
    '''Function to write SNANA fit file and label it with processor ID
    Input: nid = processor ID
    '''
    print "WRITE NML"
    Filename='ONESIM_SNFIT_DES_MASTER.NML'
    SIMLIBlines = []
    with open(Filename, 'r') as file:
        SIMLIBlines=file.read().splitlines()
    SIMLIBlines = np.array( SIMLIBlines)
    SIMLIBlines[2] = "OUTDIR:  OUT_ABC_"+str(nid)
    SIMLIBlines[4] = "VERSION: EJ_ABC_"+str(nid)
    SIMLIBlines[11] = "     MAGOBS_SHIFT_ZP = 'i "+str(round(ZPoff0,3))+" g "+str(round(ZPoff1,3))+" r "+str(round(ZPoff2,3))+" z "+str(round(ZPoff3,3))+"'"
    Outname = Filename+'_'+str(nid)
    with open(Outname, 'w') as file:
        for l in range(len(SIMLIBlines)):
            file.write("%s\n" % SIMLIBlines[l])







def launch_sim(node_id,summary_file):

    '''Function to launch SNANA sim using subprocess
    Input: node_id = processor ID
    summary_file: name out sim summary file, once this exists sim is done and we exit the func
    Returns: error message from func check_abort
    '''
    print "LAUNCHING SIM"
    cmd='sim_SNmix.pl ONESIM_MASTER.INPUT'+'_'+str(node_id)
    process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
    while os.path.isfile(summary_file) == False:
        sleep(0.1)
    print "SIM DONE"
    return check_abort(summary_file)

def check_abort(sfile):

    '''Function which checks for words "abort" or "ABORT" in sim summary file
    Input: SNANA sim summary file
    Returns: sim_ok flag = True or False
    '''

    print "CHECKING FOR 'ABORT'"
    f = open(sfile, 'r')
    contents = f.read()
    f.close()
    sim_ok = True
    searchterms = ["abort","ABORT"]
    for st in searchterms:
        searchexp = re.compile("(%s)" % st, re.M)
        if searchexp.search(contents):
          sim_ok = False
    return sim_ok

def launch_fit(node_id,sim_root_file):
    
    '''Function to launch SNANA fitting program using split_and_fit
    Input: node_id = Processor ID, sim_root_file = ROOT file from fit
    '''
    cmd='split_and_fit.pl ONESIM_SNFIT_DES_MASTER.NML_'+str(node_id)
    process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
    while os.path.isfile(sim_root_file) == False:
        sleep(0.1)
    print "FIT DONE"

def check_fit(node_id,merge_file):

    '''Function to check the output of SNANA fit completed
    Input: node_id = processor ID, merge_file = merge file
    Returns: Boolean True or False depending on if "Done." appears in merge file
    '''
    f = open(merge_file, 'r')
    contents = f.read()
    f.close()
    fit_done = False
    searchterms = ["Done."]
    for st in searchterms:
        searchexp = re.compile("(%s)" % st, re.M)
        if searchexp.search(contents):
            fit_done = True
    return fit_done

def check_order_params(param_names,theta):
    '''Checks for specific parameter names and then returns ordered list for sampler
    Input: theta - list of parameters to be varied
            param_names - parsed names for each parameter
    Returns: theta_om , theta_w0, theta_alpha, theta_beta
    '''
    if len(np.where(param_names == 'om')[0]) ==0:
        theta_om = 0.3
    else:
        theta_om = theta[np.where(param_names == "om")[0][0]]
    if len(np.where(param_names == 'w0')[0])==0:
        theta_w0 = -1.0
    else:
        theta_w0 = theta[np.where(param_names == "w0")[0][0]]
    if len(np.where(param_names == "alpha")[0])==0:
        theta_alpha = 0.14
    else:
        theta_alpha = theta[np.where(param_names == "alpha")[0][0]]
    if len(np.where(param_names == "beta")[0])==0:
        theta_beta = 3.2
    else:
        theta_beta = theta[np.where(param_names == "beta")[0][0]]
    if len(np.where(param_names == "magoff")[0])==0:
        theta_magoff = 0.0
    else:
        theta_magoff = theta[np.where(param_names == "magoff")[0][0]]
    if len(np.where(param_names == "ZPg")[0])==0:
        theta_ZPg = 0.0
    else:
        theta_ZPg = theta[np.where(param_names == "ZPg")[0][0]]
    if len(np.where(param_names == "ZPr")[0])==0:
        theta_ZPr = 0.0
    else:
        theta_ZPr = theta[np.where(param_names == "ZPr")[0][0]]
    if len(np.where(param_names == "ZPi")[0])==0:
        theta_ZPi = 0.0
    else:
        theta_ZPi = theta[np.where(param_names == "ZPi")[0][0]]
    if len(np.where(param_names == "ZPz")[0])==0:
        theta_ZPz = 0.0
    else:
        theta_ZPz = theta[np.where(param_names == "ZPz")[0][0]]
    return theta_om , theta_w0, theta_alpha, theta_beta, theta_magoff, theta_ZPg, theta_ZPr, theta_ZPi, theta_ZPz

