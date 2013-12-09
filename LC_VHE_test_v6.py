import asciitable
import numpy
import matplotlib as mpt
from matplotlib import pyplot
from VOT import *
import operator
from matplotlib.patches import Ellipse
from prettytable import *
import glob
import os
'''
Ian's NASA Summer 2013 Work
'''
'''
This is the meat and bones of the tool. 
'''
class all_bursts:
	def __init__(self, eblModel = 'Dominguez', instrument = 'VERITAS', zenith = 20):
	
		self.eblModel = eblModel
		self.instrument = instrument
		self.zenith = zenith
		GRBs = []
		myGRB = {}; self.myGRB = myGRB # establishes a dictionary of all GRBs
		# accesses list of GRBs
		directory = os.getcwd()
		os.chdir(directory + '\GRBs')
		files = list(glob.glob('*.csv'))
		
		# list of all known redshifts for the current sample of 23 bursts
		redshifts = {
		'130518580': 2.49,
		'130427324': 0.34,
		'110731465': 2.83,
		'100414097': 1.368,
		'091003191': 0.8969,
		'090926181': 2.1062,
		'090902462': 1.822,
		'090510016': 0.903,
		'090328401': 0.736,
		'090323002': 3.57,
		'080916009': 4.35}
		
		for f in files:
			GRBs.append(f[:9])
		os.chdir(directory)
		# establishes the default redshift to be 2.00 if there is no known redshift
		for g in GRBs:
			if g in redshifts:
				pass
			else:
				redshifts[g] = 2.00
			
			myGRB['{}'.format(g)] = take_data(Burst = g, eblModel = self.eblModel, instrument = self.instrument, zenith = self.zenith, redshift = redshifts[g])
			self.redshift = redshifts[g]
'''
Used to analyze individual bursts. Must have corresponding csv file for the burst in the GRB folder.
The tool requires the burst name, eblModel, redshift, instrument, and zenith angle. Default instrument is Veritas.  
'''			
class take_data:

	def __init__(self, Burst = '130427A', redshift = 0.34, eblModel = 'Dominguez', instrument = 'VERITAS', zenith = 20):	
		
		#Although it is a bit unclear to me why, you have to put a self. before things so that the code can refer to that piece of code.
		self.redshift = redshift
		self.eblModel = eblModel
		self.instrument = instrument
		self.zenith = zenith
		

	# this finds the GRB file and reads it there is another file I created that will take txt files and convert it to a more friendly csv file 
		GRB = 'GRBs/'+ Burst +'.csv'
		burst = asciitable.read(GRB, delimiter = ',')
	
	# this makes the arrays(tuple) to which we will add values from the GRB table
	
		
	
	
	# this imports the data from the csv files 
	# It looks specifically for certain header names. It's very picky so I suggest using the converter file I created.
	# Now, ideally, we have all the data we could possibly need from the csv files
		self.start_time = []
		for x_1 in burst['start time']:
			self.start_time.append(x_1)
		self.stop_time = []
		for x_2 in burst['stop time']:
			self.stop_time.append(x_2)
		self.PF = []
		for y in burst['photon flux']:
			self.PF.append(y)
		self.PF_err = []
		for y_1 in burst['photon flux err']:
			self.PF_err.append(y_1)
		self.PI = []
		for p in burst['Photon index']:
			self.PI.append(p)
		self.PI_err = []
		for p_1 in burst['Photon index error']:
			self.PI_err.append(p_1)
		self.Flux_erg = []
		for E_1 in burst['Energy flux (erg/cm2/s)']:
			self.Flux_erg.append(E_1)
		self.Flux_erg_error = []
		for E_e in burst['Energy flux error']:
			self.Flux_erg_error.append(E_e)
		
		
		
	 
	#photon index error calculator
		PIP, PIM = zip(*[(pi+pie, pi-pie) for pi, pie in zip(self.PI, self.PI_err)])
		self.PIP = PIP; self.PIM = PIM
		
	
	# these make new arrays. One is the center of each time bin and the other is the time between the center and edge of each time bin
		time,time_err = zip(*[(((y-x)/2)+x,(y-x)/2) for x,y in zip(self.start_time, self.stop_time)])
		self.time = time; self.time_err = time_err
	
	# So we have all these time bins and for loops but it would be nice, say if we wanted the differential flux for a single time bin, to call a specific time bin
	# well I was nice enough to do it for you
		time_bin = range(0, len(time), 1)

	# This changes the Flux from the LAT from ergs to GeV 
		Flux_GeV =  [624.150934 * x for x in self.Flux_erg]; self.Flux_GeV = Flux_GeV
		
		Flux_GeV_error =  [624.150934 * x for x in self.Flux_erg_error]; self.Flux_GeV_error = Flux_GeV_error
	# The effective energy range for the LAT data spans from 0.1 GeV to 100 GeV 
		E2 =  100 # GeV
		E1 = 0.1 # GeV

	# This neat little formula saved me a bit of time. Rather than do the calculation and write a new csv file, this does the calculations for the Normalization values N0 here.
		self.Norm = []; self.Norm_P = []; self.Norm_M = []
		for F, F_err, I in zip(self.PF, self.PF_err, self.PI):
			self.Norm.append(((F/(E2-E1))*((E2**(1+I)-E1**(1+I))/(1+I)))); 
			self.Norm_P.append((((F+F_err)/(E2-E1))*((E2**(1+I)-E1**(1+I))/(1+I)))); 
			self.Norm_M.append((((F-F_err)/(E2-E1))*((E2**(1+I)-E1**(1+I))/(1+I)))); 
		
	# A dictionary so that we can refer to specific time bins later
	 	myVOT = {}; self.myVOT = myVOT
		self.One_GeV = []
		self.Four_GeV = []
		self.One_TeV = []
		self.int_flux = []
		self.detTimes = []
		self.crab_flux = []
		# the nominal calculation. No errors involved
		for tb, i, v in zip(time_bin, self.Norm, self.PI):
			myVOT['{}'.format(tb)] = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			redshift = self.redshift
			
			
			# makes a boolean array (Trues and Falses) of all the energy bins True if E >  minimum safe Energy
			mymask = (myVOT['{}'.format(tb)].VS.EBins > 10**myVOT['{}'.format(tb)].VR.EASummary['minSafeE']) 
			#print np.shape(mymask),type(myVOT['{}'.format(tb)].VS.EBins),type(myVOT['{}'.format(tb)].VS.dNdE)
			#print mask
			
			# makes a new array of energies above 10**self.VR.EASummary['minSafeE'] GeV
			E_new = [10**myVOT['{}'.format(tb)].VR.EASummary['minSafeE']] + myVOT['{}'.format(tb)].VS.EBins[mymask]
			# makes a new differential flux with an interpolated value at 10**self.VR.EASummary['minSafeE'] GeV and all values above (not interpolated)
			dNdE_new = [np.interp(10**myVOT['{}'.format(tb)].VR.EASummary['minSafeE'], myVOT['{}'.format(tb)].VS.EBins, myVOT['{}'.format(tb)].VS.dNdE_absorbed)] + np.array(myVOT['{}'.format(tb)].VS.dNdE_absorbed)[mymask] 
			# interpolates the differential flux at 1 GeV, 400 GeV and 1 TeV 
			self.One_GeV.append(np.interp(1,myVOT['{}'.format(tb)].VS.EBins, myVOT['{}'.format(tb)].VS.dNdE_absorbed))
			self.Four_GeV.append(np.interp(400,myVOT['{}'.format(tb)].VS.EBins, myVOT['{}'.format(tb)].VS.dNdE_absorbed))
			self.One_TeV.append(np.interp(1000,myVOT['{}'.format(tb)].VS.EBins, myVOT['{}'.format(tb)].VS.dNdE_absorbed))
			
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			self.int_flux.append(np.trapz(dNdE_new, x = E_new)) 
			
			#checks how long detection takes
			
			crabFlux = 100*myVOT['{}'.format(tb)].rate*60./myVOT['{}'.format(tb)].VR.crabRate
			self.crab_flux.append(myVOT['{}'.format(tb)].VR.crabRate)
			detTime = np.interp([crabFlux*0.01], myVOT['{}'.format(tb)].VR.SensCurve[:,0], myVOT['{}'.format(tb)].VR.SensCurve[:,1])*60
			self.detTimes.append((detTime[0]))
		# Only taking into account the difference in positive normalization values
		myVOT_PFP = {}; self.myVOT_PFP = myVOT_PFP
		
		'''
		Calculates the upper limit of the extrapolation based on the larger normalization value. 
		This is used to create the upper portion of the flux error bars. (i.e. N0 + error)
		'''
		self.int_flux_p = []
		self.detTimes_p = []
		self.One_GeV_P = []
		self.Four_GeV_P = []
		self.One_TeV_P = []
		self.crab_flux_p = []
		for tb, i, v in zip(time_bin, self.Norm_P, self.PI):
			myVOT_PFP['{}'.format(tb)] = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energies above the minimum safe energy range
			mymask_p = (myVOT_PFP['{}'.format(tb)].VS.EBins > 10**myVOT_PFP['{}'.format(tb)].VR.EASummary['minSafeE'])
		
			
			# Limits the energies to only those within the safe energy range
			E_new_p = [10**myVOT_PFP['{}'.format(tb)].VR.EASummary['minSafeE']] + myVOT_PFP['{}'.format(tb)].VS.EBins[mymask_p]
			
			# Limits the differential flux values (i.e. y axis )
			dNdE_new_p = [np.interp(10**myVOT_PFP['{}'.format(tb)].VR.EASummary['minSafeE'], myVOT_PFP['{}'.format(tb)].VS.EBins, myVOT_PFP['{}'.format(tb)].VS.dNdE_absorbed)] + np.array(myVOT_PFP['{}'.format(tb)].VS.dNdE_absorbed)[mymask_p]
			
			# Uses the interpolation function to determine the differential flux for energies 1 GeV, 400 GeV and 1 TeV 
			self.One_GeV_P.append(np.interp(1,myVOT_PFP['{}'.format(tb)].VS.EBins, myVOT_PFP['{}'.format(tb)].VS.dNdE_absorbed))
			self.Four_GeV_P.append(np.interp(400,myVOT_PFP['{}'.format(tb)].VS.EBins, myVOT_PFP['{}'.format(tb)].VS.dNdE_absorbed))
			self.One_TeV_P.append(np.interp(1000,myVOT_PFP['{}'.format(tb)].VS.EBins, myVOT_PFP['{}'.format(tb)].VS.dNdE_absorbed))
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			self.int_flux_p.append(np.trapz(dNdE_new_p, x = E_new_p))
			
			#checks how long detection takes
			
			crabFlux_p = 100*myVOT_PFP['{}'.format(tb)].rate*60./myVOT_PFP['{}'.format(tb)].VR.crabRate
			self.crab_flux_p.append(myVOT_PFP['{}'.format(tb)].VR.crabRate)
			detTime_P = np.interp([crabFlux_p*0.01], myVOT_PFP['{}'.format(tb)].VR.SensCurve[:,0], myVOT_PFP['{}'.format(tb)].VR.SensCurve[:,1])*60
			self.detTimes_p.append(detTime_P[0])
			
			myVOT_PFM = {}; self.myVOT_PFM = myVOT_PFM
			self.int_flux_m = []
			self.detTimes_m = []
			self.One_GeV_M = []
			self.Four_GeV_M = []
			self.One_TeV_M = []
			self. crab_flux_m = []
			# only taking into account the difference in negative normalization error
		for tb, i, v in zip(time_bin, self.Norm_M, self.PI):
			myVOT_PFM['{}'.format(tb)] = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energy bins depending on whether True if E > 10**self.VR.EASummary['minSafeE'] GeV
			mymask_m =  (myVOT_PFM['{}'.format(tb)].VS.EBins > 10**myVOT_PFM['{}'.format(tb)].VR.EASummary['minSafeE'])
			#print np.shape(mymask),type(myVOT['{}'.format(tb)].VS.EBins),type(myVOT['{}'.format(tb)].VS.dNdE)
			#print mask
			
			# makes a new array of energies above the minimum sensitivity
			E_new_m = [10**myVOT_PFM['{}'.format(tb)].VR.EASummary['minSafeE']] + myVOT_PFM['{}'.format(tb)].VS.EBins[mymask_m]
			# makes a new differential flux with an interpolated value at 10**self.VR.EASummary['minSafeE'] GeV and all values above (not interpolated)
			dNdE_new_m = [np.interp(10**myVOT_PFM['{}'.format(tb)].VR.EASummary['minSafeE'], myVOT_PFM['{}'.format(tb)].VS.EBins, myVOT_PFM['{}'.format(tb)].VS.dNdE_absorbed)] + np.array(myVOT_PFM['{}'.format(tb)].VS.dNdE_absorbed)[mymask_m]		
			# interpolates the differential flux at 1 GeV, 400 GeV and 1 TeV 
			self.One_GeV_M.append(np.interp(1,myVOT_PFM['{}'.format(tb)].VS.EBins, myVOT_PFM['{}'.format(tb)].VS.dNdE_absorbed))
			self.Four_GeV_M.append(np.interp(400,myVOT_PFM['{}'.format(tb)].VS.EBins, myVOT_PFM['{}'.format(tb)].VS.dNdE_absorbed))
			self.One_TeV_M.append(np.interp(1000,myVOT_PFM['{}'.format(tb)].VS.EBins, myVOT_PFM['{}'.format(tb)].VS.dNdE_absorbed))
			
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			self.int_flux_m.append(np.trapz(dNdE_new_m, x = E_new_m))
			
			#checks how long detection takes
			
			crabFlux_m = 100*myVOT_PFM['{}'.format(tb)].rate*60./myVOT_PFM['{}'.format(tb)].VR.crabRate
			self.crab_flux_m.append(myVOT_PFM['{}'.format(tb)].VR.crabRate)
			detTime_m = np.interp([crabFlux_m*0.01], myVOT_PFM['{}'.format(tb)].VR.SensCurve[:,0], myVOT_PFM['{}'.format(tb)].VR.SensCurve[:,1])*60
			self.detTimes_m.append(detTime_m[0])
		
		self.perror = []
		for err, perr in zip(self.int_flux, self.int_flux_p):
			self.perror.append(perr - err)
		
		self.merror = []
		for err, merr in zip(self.int_flux, self.int_flux_m):
			self.merror.append(err-merr)
		
		
		'''
		Performs all calculations taking into account all possible sources of statistical error 
		(i.e. Normalization error and Photon Index error).
		'''
		myVOT_PIP = {}; self.myVOT_PIP = myVOT_PIP
		self.int_flux_pi = []
		self.detTimes_pi = []
		self.One_GeV_PI = []
		self.Four_GeV_PI = []
		self.One_TeV_PI = []
		self.crab_flux_pi = []
		#positive photon index error and positive normalization error
		for tb, i, v in zip(time_bin, self.Norm_P, self.PIP):
			myVOT_PIP['{}'.format(tb)] = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			self.index_p = v
			# makes a boolean array (Trues and Falses) of all the energy bins True if E > 10**self.VR.EASummary['minSafeE'] GeV
			mymask_pi = (myVOT_PIP['{}'.format(tb)].VS.EBins > 10**myVOT_PIP['{}'.format(tb)].VR.EASummary['minSafeE'])
		
			
			# makes a new array of energies above 10**self.VR.EASummary['minSafeE'] GeV
			E_new_pi = [10**myVOT_PIP['{}'.format(tb)].VR.EASummary['minSafeE']] + myVOT_PIP['{}'.format(tb)].VS.EBins[mymask_pi]
			# makes a new differential flux with an interpolated value at 10**self.VR.EASummary['minSafeE'] GeV and all values above (not interpolated)
			dNdE_new_pi = [np.interp(10**myVOT_PIP['{}'.format(tb)].VR.EASummary['minSafeE'], myVOT_PIP['{}'.format(tb)].VS.EBins, myVOT_PIP['{}'.format(tb)].VS.dNdE_absorbed)] + np.array(myVOT_PIP['{}'.format(tb)].VS.dNdE_absorbed)[mymask_pi]
			# interpolates the differential flux at 1 GeV, 400 GeV and 1 TeV 
			self.One_GeV_PI.append(np.interp(1,myVOT_PIP['{}'.format(tb)].VS.EBins, myVOT_PIP['{}'.format(tb)].VS.dNdE_absorbed))
			self.Four_GeV_PI.append(np.interp(400,myVOT_PIP['{}'.format(tb)].VS.EBins, myVOT_PIP['{}'.format(tb)].VS.dNdE_absorbed))
			self.One_TeV_PI.append(np.interp(1000,myVOT_PIP['{}'.format(tb)].VS.EBins, myVOT_PIP['{}'.format(tb)].VS.dNdE_absorbed))
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			self.int_flux_pi.append(np.trapz(dNdE_new_pi, x = E_new_pi))
			
			
			#checks how long detection takes
			crabFlux_pi = 100*myVOT_PIP['{}'.format(tb)].rate*60./myVOT_PIP['{}'.format(tb)].VR.crabRate
			self.crab_flux_pi.append(myVOT_PIP['{}'.format(tb)].VR.crabRate)
			detTime_PI = np.interp([crabFlux_pi*0.01], myVOT_PIP['{}'.format(tb)].VR.SensCurve[:,0], myVOT_PIP['{}'.format(tb)].VR.SensCurve[:,1])*60
			self.detTimes_pi.append(detTime_PI[0])
			
			
		myVOT_PIM = {}; self.myVOT_PIM = myVOT_PIM
		self.int_flux_mi = []
		self.detTimes_mi = []
		self.One_GeV_MI = []
		self.Four_GeV_MI = []
		self.One_TeV_MI = []
		self.crab_flux_mi = []
		# negative photon index error and normalization error
		for tb, i, v in zip(time_bin, self.Norm_M, self.PIM):
			myVOT_PIM['{}'.format(tb)] = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			self.index_m = v
			# makes a boolean array (Trues and Falses) of all the energy bins True if E > 10**self.VR.EASummary['minSafeE'] GeV
			mymask_mi =  (myVOT_PIM['{}'.format(tb)].VS.EBins > 10**myVOT_PIM['{}'.format(tb)].VR.EASummary['minSafeE'])
			#print np.shape(mymask),type(myVOT['{}'.format(tb)].VS.EBins),type(myVOT['{}'.format(tb)].VS.dNdE)
			#print mask
			
			# makes a new array of energies above 10**self.VR.EASummary['minSafeE'] GeV
			E_new_mi = [10**myVOT_PIM['{}'.format(tb)].VR.EASummary['minSafeE']] + myVOT_PIM['{}'.format(tb)].VS.EBins[mymask_m]
			# makes a new differential flux with an interpolated value at 10**self.VR.EASummary['minSafeE'] GeV and all values above (not interpolated)
			dNdE_new_mi = [np.interp(10**myVOT_PIM['{}'.format(tb)].VR.EASummary['minSafeE'], myVOT_PIM['{}'.format(tb)].VS.EBins, myVOT_PIM['{}'.format(tb)].VS.dNdE_absorbed)] + np.array(myVOT_PIM['{}'.format(tb)].VS.dNdE_absorbed)[mymask_mi]		
			# interpolates the differential flux at 1 GeV, 400 GeV and 1 TeV 
			self.One_GeV_MI.append(np.interp(1,myVOT_PIM['{}'.format(tb)].VS.EBins, myVOT_PIM['{}'.format(tb)].VS.dNdE_absorbed))
			self.Four_GeV_MI.append(np.interp(400,myVOT_PIM['{}'.format(tb)].VS.EBins, myVOT_PIM['{}'.format(tb)].VS.dNdE_absorbed)) 
			self.One_TeV_MI.append(np.interp(1000,myVOT_PIM['{}'.format(tb)].VS.EBins, myVOT_PIM['{}'.format(tb)].VS.dNdE_absorbed)) 
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			self.int_flux_mi.append(np.trapz(dNdE_new_mi, x = E_new_mi))
			
			#checks how long detection takes
			
			crabFlux_mi = 100*myVOT_PIM['{}'.format(tb)].rate*60./myVOT_PIM['{}'.format(tb)].VR.crabRate
			self.crab_flux_mi.append(myVOT_PIM['{}'.format(tb)].VR.crabRate)
			detTime_MI = np.interp([crabFlux_mi*0.01], myVOT_PIM['{}'.format(tb)].VR.SensCurve[:,0], myVOT_PIM['{}'.format(tb)].VR.SensCurve[:,1])*60
			self.detTimes_mi.append(detTime_MI[0])
	
		
		# Lots and lots of error calculations
		One_GeV_M_error,Four_GeV_M_error, One_TeV_M_error = zip(*[((x-s),(y-t),(z-u)) for s,t,u,x,y,z in zip(self.One_GeV_M, self.Four_GeV_M, self.One_TeV_M, self.One_GeV, self.Four_GeV, self.One_TeV)])
		
		self.One_GeV_M_error = One_GeV_M_error; self.Four_GeV_M_error = Four_GeV_M_error; self.One_TeV_M_error = One_TeV_M_error
		
		One_GeV_P_error,Four_GeV_P_error, One_TeV_P_error = zip(*[((s-x),(t-y),(u-z)) for s,t,u,x,y,z in zip(self.One_GeV_M, self.Four_GeV_M, self.One_TeV_M, self.One_GeV, self.Four_GeV, self.One_TeV)])
		
		self.One_GeV_P_error = One_GeV_P_error; self.Four_GeV_P_error = Four_GeV_P_error; self.One_TeV_P_error = One_TeV_P_error
		
		# percent error
	def error(self, timebin = 1):
		pp_error = (np.absolute(self.int_flux_pi[timebin] - self.int_flux[timebin])/self.int_flux[timebin])*100
		pm_error = (np.absolute(self.int_flux_mi[timebin] - self.int_flux[timebin])/self.int_flux[timebin])*100
		self.pp_error = pp_error
		self.pm_error = pm_error
		# a plot of the integrated flux
	def intplot(self):
		fig2 = pyplot.figure(figsize=(16,8))	
		fig2ax1 = fig2.add_subplot(111)
		#fig2ax2 = fig2.add_subplot(112)
		fig2ax1.set_ylabel(r'Integral Flux [cm$^{-2}$ s$^{-1}$ ]', fontsize = '20')
		fig2ax1.set_xlabel('Time [s]', fontsize = '20')
		fig2ax1.set_xscale('log')
		fig2ax1.set_yscale('log')
		msk1 = 2*np.array(self.time_err) > np.array(self.detTimes)
		msk2 = 2*np.array(self.time_err) < np.array(self.detTimes)
		l1 = fig2ax1.errorbar(np.array(self.time)[msk1], np.array(self.int_flux)[msk1], xerr = np.array(self.time_err)[msk1],
			yerr = [np.array(self.merror)[msk1], np.array(self.perror)[msk1]], fmt='_', color = 'g', label = 'Detectable')
		l11 = fig2ax1.errorbar(np.array(self.time)[msk2], np.array(self.int_flux)[msk2], xerr = np.array(self.time_err)[msk2],
			yerr = [np.array(self.merror)[msk2], np.array(self.perror)[msk2]], fmt='_', color = 'r', label = 'Detectable')
		#l2 = fig2ax1.scatter(np.array(self.time)[msk1], np.array(self.int_flux)[msk1], color = 'g')
		#l22 = fig2ax1.scatter(np.array(self.time)[msk2], np.array(self.int_flux)[msk2], color = 'r')
		fig2ax2 = fig2ax1.twinx()
		fig2ax1.legend(fontsize = '18')
		#fig2ax1.set_xlim(1,1e6)
		#fig2ax1.set_ylim(1e-12,1e-2)
		fig2ax1.axvline(100, color = 'k')
		# a plot of LAT flux, photon index, and Predicted telescope flux
	def plotcompare(self):
		fig = mpt.pyplot.gcf()
		fig.set_size_inches(18.5,10.5)
		
		ax1 = fig.add_subplot(3,1,1)
		ax2 = fig.add_subplot(3,1,2)
		ax3 = fig.add_subplot(3,1,3)
		fig.tight_layout()
		ax1.set_ylabel(r'LAT Flux [GeV$^{-1}$cm$^{-2}$ s$^{-1}$ ]', fontsize = '20')
		ax2.set_ylabel(r'Photon Index', fontsize = '20')
		ax3.set_ylabel(r'Predicted Veritas Flux [cm$^{-2}$ s$^{-1}$ ]', fontsize = '20')
		ax3.set_xlabel('Time [s]', fontsize = '30')
		#ax1.set_title('Sharing both axes')
		ax1.set_xscale('log')
		ax1.set_yscale('log')
		ax2.set_xscale('log')
		ax3.set_xscale('log')
		ax3.set_yscale('log')
		msk1 = 2*np.array(self.time_err) > np.array(self.detTimes)
		msk2 = 2*np.array(self.time_err) < np.array(self.detTimes)
		#ax3.annotate('Photon Index error ranges from',(0.8, 0.5),
                # xycoords="axes fraction", va="bottom", ha="center",
                # bbox=dict(boxstyle="round, pad=1", fc="w"))

		ax1.errorbar(self.time, self.Flux_GeV, xerr = self.time_err, yerr = self.Flux_GeV_error, fmt = '_', color = 'k')
		ax2.errorbar(self.time, self.PI, xerr = self.time_err, yerr = self.PI_err, fmt = '_', color = 'k')
		l1 = ax3.errorbar(np.array(self.time)[msk1], np.array(self.int_flux)[msk1], xerr = np.array(self.time_err)[msk1],
			yerr = [np.array(self.merror)[msk1], np.array(self.perror)[msk1]], fmt='_', color = 'g', label = 'Detectable')
		l2 = ax3.errorbar(np.array(self.time)[msk2], np.array(self.int_flux)[msk2], xerr = np.array(self.time_err)[msk2],
			yerr = [np.array(self.merror)[msk2], np.array(self.perror)[msk2]], fmt='_', color = 'r', label = 'Undetectable')
		ax3.axvline(100, color = 'k')
		el = Ellipse((2, -1), 0.5, 0.5)
		ax3.legend(fontsize = '20')
		# Fine-tune figure; make subplots close to each other and hide x ticks for
		# all but bottom plot.
		fig.subplots_adjust(hspace = 0.1)
		pyplot.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
	def bnplot(self, timebin = '5'):
		fig2 = mpt.pyplot.gcf()
		fig2.set_size_inches(18.5,10.5)
		fig2ax1 = fig2.add_subplot(111)
		fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]', fontsize = '20')
		fig2ax1.set_xlabel('Energy [GeV]', fontsize = '30')
		fig2ax1.set_ylim(1e-25, 10e6)
		fig2ax1.set_xlim(1,1e6)
		l1 = fig2ax1.loglog(self.myVOT[timebin].VS.EBins, self.myVOT[timebin].VS.dNdE, 'b', linestyle = '--')
		l2 = fig2ax1.loglog(self.myVOT_PIP[timebin].VS.EBins, self.myVOT_PIP[timebin].VS.dNdE, 'g', linestyle = '--')
		l3 = fig2ax1.loglog(self.myVOT_PIM[timebin].VS.EBins, self.myVOT_PIM[timebin].VS.dNdE, 'r', linestyle = '--')
		l1a = fig2ax1.loglog(self.myVOT[timebin].VS.EBins, self.myVOT[timebin].VS.dNdE_absorbed, 'b')
		l2a = fig2ax1.loglog(self.myVOT_PIP[timebin].VS.EBins, self.myVOT_PIP[timebin].VS.dNdE_absorbed, 'g')
		l3a = fig2ax1.loglog(self.myVOT_PIM[timebin].VS.EBins, self.myVOT_PIM[timebin].VS.dNdE_absorbed, 'r')
		
		fig2ax2 = fig2ax1.twinx()
		l0 = fig2ax2.loglog((10**self.myVOT[timebin].VR.EACurve[0:,0])[::2],(self.myVOT[timebin].VR.EACurve[0:,1])[::2],'ko')
		fig2ax2.set_ylabel(r'Effective Area [cm$^2$]', fontsize = '17')
		fig2ax2.axvline(10**self.myVOT[timebin].VR.EASummary['minSafeE'],color='k')
		fig2ax2.axvline(10**self.myVOT[timebin].VR.EASummary['maxSafeE'],color='k')
		fig2.legend([l1a[0], l2a[0], l3a[0], l1[0], l2[0], l3[0], l0[0]], ['Nominal value','Positive photon index error','Negative photon index error','Nominal unabsorbed','Positive unabsorbed', 'negative unabsorbed', 'EA'], 1, fontsize = '18')
		pyplot.show()
	def dnplot(self):
		fig2 = pyplot.figure(figsize=(16,8))
		fig2ax1 = fig2.add_subplot(111)
		fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
		fig2ax1.set_xlabel('Time [s]')
		fig2ax1.set_xscale('log')
		fig2ax1.set_yscale('log')
		timer = np.arange(1, 10000, 1)
		x = 1 # 1 GeV
		y = 400 # 400 GeV
		z = 1000 # 1000 GeV
		crab_flux_1 = 2.83e-11*(x/1000.)**(-2.62) 
		crab_flux_2 = 2.83e-11*(y/1000.)**(-2.62) 
		crab_flux_3 = 2.83e-11*(z/1000.)**(-2.62) 
		cflux1 = []
		cflux2 = []
		cflux3 = []
		for t in timer:
   			cflux1.append(crab_flux_1)
    		cflux2.append(crab_flux_2)
    		cflux3.append(crab_flux_3)
		
		l20 = fig2ax1.errorbar(np.array(self.time), np.array(self.One_GeV), xerr = self.time_err, yerr = [np.array(self.One_GeV_M_error), np.array(self.One_GeV_P_error)], fmt = '+', color = 'b')
		l21 = fig2ax1.errorbar(np.array(self.time), np.array(self.Four_GeV), xerr = self.time_err, yerr = [np.array(self.Four_GeV_M_error), np.array(self.Four_GeV_P_error)], fmt = '+', color = 'g')
		l22 = fig2ax1.errorbar(np.array(self.time), np.array(self.One_TeV), xerr = self.time_err, yerr = [np.array(self.One_TeV_M_error), np.array(self.One_TeV_P_error)], fmt = '+', color = 'r')
		l31 = fig2ax1.plot(timer, cflux1)
		l32 = fig2ax1.plot(timer, cflux2)
		l33 = fig2ax1.plot(timer, cflux3)
		
	
	
		fig2.legend([l20[0],l21[0],l22[0]], ['1 GeV', '400 GeV','1 TeV '], 1)
		
		fig2ax1.set_xlim(1e-2,1e6)
		fig2ax1.set_ylim(1e-15,1e4)
		fig2ax1.axvline(100, color = 'k')
		pyplot.show()
			
	def make_table(self, timebin = 6):
		if 2*self.time_err[timebin] >= self.detTimes[timebin]:
			detection = 'Detectable'
		else:
			detection = 'Not Detectable'
		if 2*self.time_err[timebin] >= self.detTimes_pi[timebin]:
			detection_p = 'Detectable'
		else:
			detection_p = 'Not Detectable'
		if 2*self.time_err[timebin] >= self.detTimes_mi[timebin]:
			detection_m = 'Detectable'
		else:
			detection_m = 'Not Detectable'
		self.detection = detection	
		self.detection_p = detection_p
		self.detection_m = detection_m
		pp_error = (np.absolute(self.int_flux_pi[timebin] - self.int_flux[timebin])/self.int_flux[timebin])*100
		pm_error = (np.absolute(self.int_flux_mi[timebin] - self.int_flux[timebin])/self.int_flux[timebin])*100
		self.pp_error = pp_error
		self.pm_error = pm_error
		t = PrettyTable(["Photon Index","Integrated Flux","Percent Error", "detection"])
		#t.align["Photon Index"] = "l" #left align photon indicies  
		t.add_row([self.PI[timebin], self.int_flux[timebin], 'N/A', self.detection])
		t.add_row([self.PIP[timebin], self.int_flux_pi[timebin], str(self.pp_error)+'%', self.detection_p])
		t.add_row([self.PIM[timebin], self.int_flux_mi[timebin], str(self.pm_error)+'%', self.detection_m])
		print t
	
	
		