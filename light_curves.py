import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
import sys
from time import sleep
from distance_calc import *
import rootpy
import re
from rootpy.io import root_open, File
from time import sleep

#This module contains the LightCurves class for reading and cleaning output from ROOT
#and functions lookup_z_mjd(), get_obs_lc(),get_smoothed_lc() and create_spline()

__all__ = ['LightCurves']

BandNames=['g','r','i','z']

class LightCurves(object):
	def __init__(self,name,fromfiles=True,ROOTfile=None,sdata=None,data=None,cosmo_params=(0.3,0.7,0.7,-1.,0.),dump_file=None,simlib=None,wmode=-1,beta=None,alph=None):
                self.name=name
                self.beta_param=beta
                self.alpha_param=alph
		if fromfiles:
			self.ROOTfile=ROOTfile
                        self.simlib = simlib
                        if self.simlib:
                            self.read_simlib()
			self.read_root()
                        if dump_file:
                            self.read_dump(dump_file)
			self.name = name
		else:
			self.data=data
			self.sdata=sdata
                self.wmode = wmode
		self.om,self.ol,self.h0,self.w0,self.wa = cosmo_params
		self.dist()

        def read_simlib(self):
            SIMLIBlines = []
            with open(self.simlib, 'r') as file:
                SIMLIBlines=file.read().splitlines()
                SIMLIBlines = np.array(SIMLIBlines)
            self.simlib_cids = [];self.simlib_libid = []
            for ii,l in enumerate(SIMLIBlines):
                if l[0] == "L" and l[1] == "I":
                    self.simlib_libid.append(int(l.split(':')[1]))
                elif l[0] == 'O' and l[1] == 'R':
                    self.simlib_cids.append(int(l.split(':')[2]))


	def dist(self):
		DC = DistanceCalc(self.om,1.-(self.om+self.ol),self.ol,self.wmode,(self.w0,self.wa),self.h0)
		self.dl_z = np.arange(0.,1.,0.001)
		self.dl	=np.array(map(DC.d_l,self.dl_z))
		self.dl_spline = UnivariateSpline(self.dl_z,self.dl,k=3)
                self.mu_func =DC.mu

	def lookup_z_mjd(self,cid):
		'''Quick lookup for z and pkmjd for one SN'''
		return self.sumdata['simZ'][self.sumdata['CIDS'] == cid], self.sumdata['simPKMJD'][self.sumdata['CIDS'] == cid]

	def lookup_z(self,cid):
		'''Quick lookup for z for one SN'''
		return self.sumdata['simZ'][self.sumdata['CIDS'] == int(cid)]


	def get_obs_lc(self,SNcid,dataflag):
		'''
		given cid for 1 SN and bandname (str), returns array where DATAFLAG == 1 with columns # MJD  Flux flux_err 
		'''
                vals= self.lcs[SNcid]
                gband = np.where((vals[5] == 2) & (vals[4] == dataflag) & (vals[6] == 0))[0]
                rband = np.where((vals[5] == 3) & (vals[4] == dataflag) & (vals[6] == 0))[0]
                iband = np.where((vals[5] == 4) & (vals[4] == dataflag) & (vals[6] == 0))[0]
                zband = np.where((vals[5] == 5) & (vals[4] == dataflag) & (vals[6] == 0))[0]
                pkmjd = vals[0][0]
                mjd = [pkmjd+vals[1][gband],pkmjd+vals[1][rband],pkmjd+vals[1][iband],pkmjd+vals[1][zband] ]
                F = [vals[2][gband],vals[2][rband],vals[2][iband],vals[2][zband] ]
                Ferr = [vals[3][gband],vals[3][rband],vals[3][iband],vals[3][zband] ]
		return mjd,F,Ferr


        def read_root(self):
            print "  ", self.ROOTfile
            waiting_printout = True 
            f=File(self.ROOTfile,"read")
            fitres_tree = None
            fitres_tree =  f.Get('FITRES')
            lc_tree = f.Get('SNLCPAK')
            self.fitres = {}
            self.NumSN = fitres_tree.GetEntries()
            if self.name =='data':
                for j ,event in enumerate(fitres_tree):
                    self.fitres[event.CID] = [event.SIM_PKMJD, event.SIM_ZCMB,j,event.x1,event.x1ERR,event.c,event.cERR,event.mB,event.mBERR]
            else:
                for j ,event in enumerate(fitres_tree):
                    #self.fitres[event.CID] = [event.SIM_PKMJD, event.SIMZ,event.SIM_LIBID,j]
                    self.fitres[event.CID] = [event.SIM_PKMJD, event.SIM_ZCMB,event.SIM_LIBID,j,event.SIM_x1,event.SIM_c,event.SIM_DLMAG,event.SIM_mB, event.c,event.cERR,event.x1,event.x1ERR,event.mB,event.mBERR]
            #read lcs
            self.lcs={}
            for i,ll in enumerate(lc_tree):
                temp =[np.array(ll.BAND_PEAKMJD),np.array(ll.TOBS),np.array( ll.FLUXCAL),np.array(ll.FLUXCAL_ERR),np.array(ll.IFLAGDATA),np.array(ll.IFILTOBS), np.array(ll.REJECT)]
                self.lcs[i] = temp


	
	def read_lcplot(self):
		dtype=[('CIDS',np.int32),('MJD',np.float32),('Tobs',np.float32),('FLUXCAL',np.float32),('FLUXCAL_ERR',np.float32),('DATAFLAG',np.int32),('BAND','|S15')]
        	self.data=np.loadtxt(self.LCPLOTfile,skiprows = 2, comments=("NVAR","VARNAMES"),dtype =dtype, usecols=(1,2,3,4,5,6,7))
		
		
	def read_dump(self,dfile):
		dtype=[('CIDS',np.int32),('libid',np.int32), ('ngen',np.int32)]
        	self.dumpdata=np.loadtxt(dfile,skiprows = 2 ,usecols = (1,6,12),dtype=dtype,comments=("NVAR","VARNAMES"))
		

	def create_spline(self,x=None,y=None,z=None):
		'''creates spline to full LC given cid for SN and band
		   input: 
			x,y,z array for spline (optional) if not provided spline is for MJD and FLUXCAL
		   returns:
			spline for LC   
		'''
                order = np.argsort(x)
                self.LCspline = interp1d(x[order],y[order])	
		return self.LCspline
