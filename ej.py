"""
Total ejecta mass calculation for the slowly-declining 01ay.
To vet the code against what Max used for his analysis


Everything presumes a 19d rise time
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from pack import fid_time,bol 
from scipy.interpolate import interp1d

def peak(lc):

    l=np.linspace(lc[:,0].min(), lc[:,0].max(), 100)
    spl=interp1d(lc[:,0], lc[:,1], kind='cubic')
    gl = spl(l)
    return max(gl), l[gl==max(gl)][0]

def chisq_by_hand(f, t0, obs):
	"""
	chi2 = (model - flux)**2/err**2
	"""
	edep = f(obs[:,0], t0)
	
	flux = obs[:,1]
	
	chisq = np.sum((edep-flux)**2/obs[:,2]**2)
	return chisq

def fit_chisq(f, obs, tarr=np.linspace(0, 50, 100)):
	chiarr = np.array([chisq_by_hand(f, i, obs) for i in tarr])
	return tarr[chiarr==min(chiarr)][0], min(chiarr)

#doesnt require a function form yet
bp = bol.bol_func().bolpeak
bolometric_lightcurve_file = sys.argv[1]#'stritzinger/sn01ay.edep.new.dat'

bollc = np.loadtxt(bolometric_lightcurve_file)
#bollc[:,1]=pow(10, bollc[:,1])
mmax, tmax = peak(bollc)

#crude if loop to not shift light curve if already shifted to phase
if max(bollc[:,0]) < 500:
    tmax = 0
    mmax = max(bollc[:,1])

bollc[:,0]-=tmax
print "The bolometric peak is:", tmax
tail = bollc[(bollc[:,0] > 40) & (bollc[:,0] < 100)]

risetime= float(sys.argv[2])
coef = bol.arn_coef(risetime)
mni = mmax/(coef*1e43)
ft = fid_time.fid_time(mni)

print "The inferred Nickel mass is:", mni

popt, pcov = curve_fit(ft.edp_nomc, tail[:,0], tail[:,1], sigma=tail[:,2], p0=[20.])

t0= popt[0]+risetime
print "fiducial timescale, t0 is:", t0, pcov[0]
print "Ejecta mass is ", ft.ejm(t0), ft.ejm_mc((t0, pcov[0]+3))

t=np.linspace(0, 200, 500)
print fit_chisq(ft.edp_nomc, tail), chisq_by_hand(ft.edp_nomc, 9, tail)



#plotting the results
plt.figure(1)
plt.errorbar(bollc[:,0]+risetime, bollc[:,1],bollc[:,2], fmt='rs')
plt.plot(t, ft.edp_nomc(t, t0), label='t0='+str(t0))
plt.plot(t, ft.edp_nomc(t, 38), label='t0=28')

plt.legend(loc=0)
plt.xlabel('Phase (days)')
plt.ylabel('Bolometric Flux')

#chisq plot
plt.figure(2); tarr=np.linspace(0, 50, 100)
plt.plot(tarr, np.array([chisq_by_hand(ft.edp_nomc,i, tail) for i in tarr]))
plt.show()
