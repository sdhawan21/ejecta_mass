"""
Code to calculate Ejecta masses of SNIa
.. - Fit an energy deposition curve to data between +40 and +100 d
.. - Use equation 24 and equation A16 of Jeffrey 1999
.. - Use values for e-folding velocity, absorption opacity and q-factor (traces Ni distribution)

This code follows the procedure in Stritzinger et al. 2006 and therefore would expectedly given out larger errors than those in Scalzo et al. 2014 )

Error from Reddening and Distance is propagated in the bolometric light curve (both for the tail and the peak that gives out the 56Ni mass)
Propagated using a Monte Carlo

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
    Note: The input Ni mass can be a realisation from the normal distribution from the peak and the error on the peak 

    """

    def __init__(self, mni):
        """
        Initiate the class with a Nickel mass measurement
        """
        self.mni=mni
        #self.e_mni = e_mni
    
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

    #list of SNe files to be extract
    filelist=np.loadtxt(sys.argv[1], dtype='string')

    #pre: is the directory with the light curves 
    #suf: is the suffix for each file (u->H coverage)

    pre='lcbol_distrib/'
    suf='_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
    arr=[]
    
    #loop over the SNe in the input file
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
	   

	    errfunc = bol.bol_func().err_peak
            
            #e_mni = np.std([errfunc(bollc) for k in range(10000)])
            e_mni = bollc[bollc[:,1]==max(bollc[:,1])][0][2]/2e43

            #very, very crude evaluation of Dm15 
	    dm15 = - (mni -1.3) /0.62
	    
            #check whether there is pre-maximum bolometric coverage
            assert tmax != bollc[0][0]
            
            #shift time axis to phase
            bollc[:,0]-=tmax
            
            #define the phase range for the curve fit (+40 to +100)
            b=bollc[(bollc[:,0]>40.) & (bollc[:,0] < 100.)]
    
            print 'Curve fitting', i[0]
            
            #best curvefit parameters
            param_arr=[]
            
            
            k=0
            # run the Monte Carlo for each realisation of 56Ni mass
            while k < 1000:
            	real = np.random.normal(mni, e_mni)
            	if real > 0 and real < 2:
            		#initialize the class with 56Ni mass realisation
            		fid=fid_time(real)

                	#optimal fitting parameters and covariance matrix 
            		popt, pcov = curve_fit(fid.edp_nomc, b[:,0], b[:,1])
          
                	#
	        	#perr = np.sqrt(np.diag(pcov))
            		if popt[0] > 0 and popt[0] < 40:
				param_arr.append(popt[0])
           			k+=1
           		
            #calculate rise time using ganeshalingham et al. 
            rt=16.5-5.*(dm15-1.1)
        
            #get best fit value (for the Mni)
            fid=fid_time(mni)
            p,c=curve_fit(fid.edp_nomc, b[:,0], b[:,1])

	    #error array from covariance matrix
            perr = np.sqrt(np.diag(c))

            # use the Monte Carlo ejecta mass function to get error from the fit
            err_fit = fid.ejm_mc([p[0], perr[0]])[1]


            #transparency time (from explosion)
            t0=p[0]+rt
        
            #error in ejecta mass from curve fit

	    ejarr = [fid.ejm(j+rt) for j in param_arr]
            ej, e_ej = fid.ejm(t0), np.std(ejarr)
	    
       	    print ej, e_ej   

            #add the error from Nickel mass and t0 fit in quadrature
            tot_err = np.sqrt(e_ej**2 + err_fit**2)

            #print "For", i[0], "the error in the peak is", perr[0]
            
            #five column array of name,  nickel mass, t0(post max), t0, mej, err    
            arr.append([i[0], round(mni, 3), round(e_mni, 3), np.mean(param_arr), t0, round(ej, 3), round(tot_err, 3)])
	
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


