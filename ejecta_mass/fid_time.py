import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.optimize import curve_fit
from pack import bol

ds2=86400.0			#number of seconds in a day
lni=1/(8.8)
lco=1/(111.3)	
sc=1e10	
fac=624150.647996		#Mev to ergs
class fid_time:
    def __init__(self, mni):
        self.mni=mni
    def edp(self, t, t0):
        Nni=self.numb(self.mni)
        t1=lni*Nni*np.exp(-lni*t)*(1.75)		#ENi
	t2=lco*Nni*(lni/(lni-lco))*(np.exp(-lco*t)-np.exp(-lni*t))	#ECo_e+
	t3=(0.12)+(3.61)*(1-np.exp(-(t0/t)**2))		#ECo_gamma
	return (t1+t2*t3)/(fac*ds2)#(t2*t3)/(fac*ds2)
    def numb(self, mni):
        return mni*1.98e33*6.022e23/56.0
    def ejm(self, t0):
        """
        Using Equation 4 of scalzo et al. to convert t0 into ejecta mass
        """
        return 8*np.pi*(((t0)*ds2)**2)*9e6*sc*3/(0.025*1.98e33)
def main():
    filelist=np.loadtxt(sys.argv[1], dtype='string')
    pre='lcbol_distrib/'
    suf='_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
    arr=[]
    for i in filelist:
        try:
            f=pre+i[0]+suf
            bollc=np.loadtxt(f, usecols=(0, 1))
            #mni=float(i[1])
            bp=bol.bol_func().bolpeak
            tmax=bp(f)[1]

            #derive Ni mass from the light curve (Arnett's rule)
            mni=bp(f)[0]/2e43

            #check whether there is pre-maximum bolometric coverage
            assert tmax != bollc[0][0]
            bollc[:,0]-=tmax
            fid=fid_time(mni)
            b=bollc[(bollc[:,0]>50.) & (bollc[:,0] < 100.)]
            popt, pcov = curve_fit(fid.edp, b[:,0], b[:,1], p0=[10.])
            arr.append([i[0], mni, popt[0], fid.ejm(popt[0]+16.5)])

            #plot the light curve and the fit to late data
            plt.plot(bollc[:,0], bollc[:,1], 'r.')
            plt.plot(b[:,0], fid.edp(b[:,0], popt[0]))
            plt.show()
        except:
            i
    np.savetxt('ejecmass.txt', arr, fmt='%s')
if __name__=="__main__":
    main()

