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
from scipy.interp1d import interpolate

def peak(lc):

    l=np.linspace(lc[:,0].min(), lc[:,0].max(), 100)
    spl=interp1d(lc[:,0], lc[:,1], kind='cubic')
    gl = spl(l)
    return max(gl), l[gl==max(gl)]

#doesnt require a function form yet
bp = bol.bol_func().bolpeak
bolometric_lightcurve_file = sys.argv[1]#'stritzinger/sn01ay.edep.new.dat'

bollc = np.loadtxt(bolometric_lightcurve_file)

mmax, tmax = peak(bollc)
bollc[:,0]-=tmax

tail = bollc[(bollc[:,0] > 50) & (bollc[:,0] < 100)]

ft = fid_time.fid_time(mmax/2e43)

popt, pcov = curve_fit(ft.edp_nomc, tail[:,0], tail[:,1])

t0= popt[0]+19
print t0, ft.ejm(t0), popt
