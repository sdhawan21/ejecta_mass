import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.integrate import quad
from scipy.optimize import curve_fit 

from pack import bol
def co_decay(t, tc):
    a = 1 - 0.967 * np.exp(-0.0347*(tc/t)**2)
    return a

def load_lc(sn, path = '/Users/lapguest/newbol/bol_ni_ej/lcbol_distrib/'):
    fil=path + sn + '_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
    lc = np.loadtxt(fil)
    tm = bol.bol_func().bolpeak(fil)[1]
    
    lc[:,0]-=tm
    return lc

def main():
    sn = sys.argv[1]
    lc=load_lc(sn)
    lc1 = lc[(lc[:,0] > 50)]
    popt, pcov = curve_fit(co_decay, lc[:,0], lc[:,1]/1e43 )

    print popt
    

main()
