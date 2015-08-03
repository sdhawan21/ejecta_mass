"""
Code to calculate Ejecta masses of SNIa
.. - Fit an energy deposition curve to data between +40 and +100 d
.. - Use equation 24 and equation A16 of Jeffrey 1999
.. - Use values for e-folding velocity, absorption opacity and q-factor (traces Ni distribution)


Note: requires pack to be install on the machine 
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.optimize import curve_fit
from scipy.integrate import quad
from pack import bol

"""
Define conversion factors and half-lives (e-folding times) of Co and Ni

"""


ds2=86400.0			#number of seconds in a day
lni=1/(8.8)
lco=1/(111.3)	
sc=1e10	
fac=624150.647996		#Mev to ergs
class fid_time:
    """
	calculate the fiducial time for a given bolometric light curve (Ni mass is obtained using Arnett's rule)
    """
    def __init__(self, mni):
        """
        Initiate the class with a Nickel mass measurement
        """
        self.mni=mni
        #self.e_mni = e_mni
    
    def edp_mc(self, t, t0):    
    	"""
    	Do Monte Carlo over the space of Mni values 
    	
    	Energy deposition function. eqn 2 Stritzinger 2006 (also see Nadyozhin 1994 for a more detailed description of the derivation)
    	
    	"""
    	
    	#if errors need to be calculated by Monte Carlo
    	
	
	Nni=self.numb(self.mni)
	
      
        t1=lni*Nni*np.exp(-lni*t)*(1.75)		#ENi
	
	t2=lco*Nni*(lni/(lni-lco))*(np.exp(-lco*t)-np.exp(-lni*t))	#ECo_e+
	t3=(0.12)+(3.61)*(1-np.exp(-(t0/t)**2))		#ECo_gamma

	return (t1+t2*t3)/(fac*ds2)#(t2*t3)/(fac*ds2)
	
	
    
    def edp_nomc(self, t, t0):
    	"""
    	Energy deposition function. eqn 2 Stritzinger 2006 (also see Nadyozhin 1994 for a more detailed description of the derivation)
    	
    	"""
    	#if errors need to be calculated by Monte Carlo
    	
	Nni = self.numb(self.mni)
	
        t1=lni*Nni*np.exp(-lni*t)*(1.75)		#ENi
	t2=lco*Nni*(lni/(lni-lco))*(np.exp(-lco*t)-np.exp(-lni*t))	#ECo_e+
	t3=(0.12)+(3.61)*(1-np.exp(-(t0/t)**2))		#ECo_gamma
	return (t1+t2*t3)/(fac*ds2)

	
	
    def numb(self, mni):
    	"""
    	Convert Ni mass to number of Ni atoms
    	"""
        return mni*1.98e33*6.022e23/56.0
    
    
    
    def ejm(self, t0):
        """
        Using Equation 4 of scalzo et al. to convert t0 into ejecta mass
        
        (talk to Wolfgang H. about this equation)
        """
        
        return 8*np.pi*(((t0)*ds2)**2)*9e6*sc*3/(0.025*1.98e33)
        
        
    def ejm_mc(self, t, n=10000):
    
    	real = np.random.normal
    
    	ar=[self.ejm(real(t[0], t[1])) for k in range(n)]
    	return np.mean(ar), np.std(ar) 
    
    def val_calc(self, sn, p):
    
	"""
	A standalone function to calculate the fiducial time and ejecta mass for a given bolometric light curve
	"""
    
    	filename = p+'bol_ni_ej/lcbol_distrib/'+sn+'_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
    	
    	t=bol.bol_func().bolpeak(filename)
    	bollc=np.loadtxt(filename)
    	bollc[:,0]-=t[1]
    	b=bollc[(bollc[:,0]>40.) & (bollc[:,0] < 100.)]
    	mni=t[0]/2e43
    	popt, pcov = curve_fit(self.edp, b[:,0], b[:,1])
    	dm15= -(mni - 1.3)/.62
    	rt=16.5-5.*(dm15-1.1)
    	return popt[0], self.ejm(popt[0]+rt)
    	
class q_fac:
	def __init__(self, q):
		self.q = q

	def func(self, z):
		return z*np.exp(-z)
		
	def denom(self, ran=[0, np.inf]):	
		d=quad(self.func, ran[0], ran[1])
		return self.q/d[0]

   
def main(path):

    filelist=np.loadtxt(sys.argv[1], dtype='string')

    #ipath='/Users/lapguest/newbol/bol_ni_ej/'
    
    #define the file prefix and suffix

    #path='/Users/lapguest/newbol/bol_ni_ej/'

    pre='lcbol_distrib/'
    suf='_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
    arr=[]
    for i in filelist:
        try:
            #complete filename 
            f=path+pre+i[0]+suf
    
    	    #light curve 
            bollc=np.loadtxt(f, usecols=(0, 1, 2))
            #mni=float(i[1])\
            
            #define functions for the peak and error of the bolometric lc
            
            bp=bol.bol_func().bolpeak
            
            err_func= bol.bol_func().err_peak
            
            tmax=bp(f)[1]
	    
	   
            #derive Ni mass from the light curve (Arnett's rule)
            
            mni=bp(f)[0]/2e43		#needs to be a function of the rise time
	   
	    e_mni = bp(f)[2]/2e43
	    
	    dm15 = - (mni -1.3) /0.62
	    
            #check whether there is pre-maximum bolometric coverage
            assert tmax != bollc[0][0]
            
            bollc[:,0]-=tmax
            
            print e_mni
            #define the phase range for the curve fit
            b=bollc[(bollc[:,0]>40.) & (bollc[:,0] < 100.)]
            #print max(b[:,0])
            print 'Curve fitting', i[0]
            
            #best curvefit parameters
            param_arr=[]
            
            for k in range(100000):
            	
            	fid=fid_time(np.random.normal(mni, e_mni))
            	popt, pcov = curve_fit(fid.edp_nomc, b[:,0], b[:,1], p0=[20.])
            #    print popt
          
	        perr = np.sqrt(np.diag(pcov))
            
		param_arr.append(popt[0])
           
            #calculate rise time using ganeshalingham et al. 
            rt=16.5-5.*(dm15-1.1)
        
            fid=fid_time(np.random.normal(mni, e_mni))
            t0=np.mean(param_arr)+rt
        
            #error in ejecta mass from curve fit

	    ejarr = [fid.ejm(j+rt) for j in param_arr]
            ej, e_ej = np.mean(ejarr), np.std(ejarr)#fid.ejm_mc([t0, perr[0]])
	           
	    print ej, e_ej   
            #print "For", i[0], "the error in the peak is", perr[0]
            
            #five column array of name,  nickel mass, t0(post max), t0, mej, err    
            arr.append([i[0], mni, np.mean(param_arr), t0, round(ej, 3), round(e_ej, 3)])
	
            #plot the light curve and the fit to late data
            """
            plt.plot(bollc[:,0], bollc[:,1], 'r.')
            #plt.plot(b[:,0], fid.edp(b[:,0], popt[0]))
            #plt.show()
	    """
        except:
            i
    np.savetxt(sys.argv[2], arr, fmt='%s')

if __name__=="__main__":
	
	if len(sys.argv) == 4:
    		
    		main(sys.argv[3])
    	else:
    		print "Usage: python", sys.argv[0], '<input file> <outfilename> <path to lcs>'
		print "Example: python", sys.argv[0], '<fileinp_ejecmass.txt> <ejecmass> </home/sdhawan>'


