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

tail = bollc[(bollc[:,0] > 40) & (bollc[:,0] < 100)]

risetime= float(sys.argv[2])
coef = bol.arn_coef(risetime)
mni = mmax/(coef*1e43)
ft = fid_time.fid_time(mni)

print "The inferred Nickel mass is:", mni

popt, pcov = curve_fit(ft.edp_nomc, tail[:,0], tail[:,1], sigma=tail[:,2])

t0= popt[0]+risetime
print "fiducial timescale, t0 is:", t0, pcov[0]
print "Ejecta mass is ", ft.ejm(t0), ft.ejm_mc((t0, pcov[0]+3))

t=np.linspace(0, 200, 500)
plt.errorbar(bollc[:,0]+risetime, bollc[:,1],bollc[:,2], fmt='rs')
plt.plot(t, ft.edp_nomc(t, t0))
plt.show()
