import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.integrate import quad
from scipy.optimize import curve_fit 

from pack import bol
def co_decay(t, tc):
    """
    Co decay function for late times (after mulytiple e-folding of 56Ni)
    Function from J. Sollerman, also used in Sollerman et al. 2004 (for 2000cx)
    """
    a =  1 - 0.967 * np.exp(-0.0347*(tc/t)**2)
    return a
def nimass_co_decay(t, t1, mni):
    a =  1 - 0.967 * np.exp((t1/t)**2)
    b = np.exp(-t/111.3)
    return 1.3e43*mni*b*a

def load_lc(sn, path = '/Users/lapguest/newbol/bol_ni_ej/lcbol_distrib/'):
    fil=path + sn + '_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
    lc = np.loadtxt(fil)
    tm = bol.bol_func().bolpeak(fil)[1]
    
    lc[:,0]-=tm
    return lc

def main():
    sn = sys.argv[2]
    lc=load_lc(sn)
    lc1 = lc[(lc[:,0] > 50)]
    popt, pcov = curve_fit(co_decay, lc[:,0], lc[:,1]/1e43 )
    
    pni, covni = curve_fit(nimass_co_decay, lc[:,0], lc[:,1])
    print pni
    fout=open('sn_tc.lis', 'a')
    fwrite = sys.argv[1]
    if fwrite:
        fout.write(sn+'\t'+str(popt)+'\n')
    fout.close()
    print popt

main()
