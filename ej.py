"""
A self-contained code to calculate the transparency timescale (and hence the ejecta mass) from the bolometric light curve tail (e.g. Stritzinger et al. 2006, Scalzo et al. 2014)

Total ejecta mass calculation for the slowly-declining 01ay. Done on bolometric light curve used in Krisciunas et al. 2011.
To vet the code against what Max used for his analysis

17-02: Includes the measured v-si values

Rise time can be entered from the command line
Usage: python ej.py <bolometric light curve filename> <rise time (days)> <initial phase(days)><plot figure or not>
Example python ej.py temp.dat 18 40 noplot
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import bol


from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
#from pack import fid_time,bol 
#from pack import bol
from sublum import fid_time

def peak(lc):
    """
    simple spline fitter for bolometric light curves 

    Args: three column array

    """

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

def fit_chisq(f, obs, tarr=np.linspace(0, 50, 1000)):
	chiarr = np.array([chisq_by_hand(f, i, obs) for i in tarr])


def chisq_by_hand(f, t0, obs, err=False):
	"""
	chi2 = (model - flux)**2/err**2
	"""
	
        edep = f(obs[:,0], t0)
	
	flux = obs[:,1]
	
	if err:
            chisq = np.sum((edep-flux)**2/obs[:,2]**2)
        else:
            chisq = np.sum((edep-flux)**2)
           
	return chisq

def fit_chisq(f, obs, tarr=np.linspace(0, 50, 1000)):
	chiarr = np.array([chisq_by_hand(f, i, obs, err=True) for i in tarr])
	return tarr[chiarr==min(chiarr)][0], min(chiarr)

#doesnt require a function form yet
bp = bol.bol_func().bolpeak
bolometric_lightcurve_file = sys.argv[1]#'stritzinger/sn01ay.edep.new.dat'

bollc = np.loadtxt(bolometric_lightcurve_file)
#bollc[:,1]=pow(10, bollc[:,1])
mmax, tmax = peak(bollc)

#crude if condition to not shift light curve if already shifted to phase
if max(bollc[:,0]) < 500:
    tmax = 0
    mmax = max(bollc[:,1])

#velocity file
vel = np.loadtxt('csptable/table_vel_max.tex', dtype='string', delimiter='&', skiprows=1, usecols=(0, -3))
vel_tab = np.array([['SN'+i[0].strip(), float(i[1][2:8])] for i in vel])
#note, no errors on the v_si

top = bolometric_lightcurve_file.find('SN')
sn_name = bolometric_lightcurve_file[top:-46]

#print vel_tab
#sys.exit()
#shift to maximum light
bollc[:,0]-=tmax

print "The bolometric peak is:", tmax

#define the late light curve as between +40 - +90 d 
p1 = int(sys.argv[3])
tail = bollc[(bollc[:,0] >= p1) & (bollc[:,0] <= 90)]
print "Number of observations in the tail is:", len(tail)
rt=True
#risetime from command line
risetime= float(sys.argv[2])
coef = bol.arn_coef(risetime)
if rt:
    mni = mmax/(coef*1e43)
else:
    mni = mmax/2e43

#
ft = fid_time.fid_time(mni)
print "The inferred Nickel mass is:", mni


popt, pcov = curve_fit(ft.edp_nomc, tail[:,0], tail[:,1], sigma=tail[:,2], p0=[20.])

#
popt, pcov = curve_fit(ft.edp_nomc, tail[:,0]+risetime, tail[:,1], sigma=tail[:,2], p0=[20.])


#the output t0 is wrt reference epoch, hence, not adding rise time
t0= popt[0]#+risetime

#calculate the v_e as V_si/sqrt(12)
#vsi = float(vel_tab[vel_tab[:,0] == sn_name][0][1])
#vecalc = vsi/np.sqrt(12)
#sys.exit()

#print "The vsi and ve_value is", vsi, vecalc


#print the fit outputs from curve_fit
print "fiducial timescale, t0 is:", t0, pcov[0]
print "Ejecta mass is ", ft.ejm(t0, calcerr=False) #ft.ejm_mc((t0, pcov[0]+3))

#define a template time axis to get the fit deposition curve
t=np.linspace(0, 200, 500)

print fit_chisq(ft.edp_nomc, tail), chisq_by_hand(ft.edp_nomc, 9, tail)


#evaluate the chisq by hand
print "The t0 and chisq minimum"
tail[:,0]+=risetime
tfit, chimin = fit_chisq(ft.edp_nomc, tail) #, chisq_by_hand(ft.edp_nomc, 28, tail)
print tfit, chimin, chimin/len(tail[:,0])


#stops at the calculation of Ejecta mass, prints velocity as well.
if sys.argv[4] !="plot":
	sys.exit()
#plotting the results
plt.figure(1)

plt.errorbar(bollc[:,0]+risetime, bollc[:,1],bollc[:,2], fmt='rs')

#plot the curves for t0, complete escape and trapping
plt.plot(t, ft.edp_nomc(t, t0)/1e43, label='t0='+str(round(t0,2)))
plt.plot(t, ft.edp_nomc(t, 0)/1e43, label='Complete Escape')
plt.plot(t, ft.edp_nomc(t, 100000)/1e43, label="Complete trapping")

plt.errorbar(bollc[:,0]+risetime, bollc[:,1]/1e43,bollc[:,2]/1e43, fmt='rs')
#plt.errorbar(tail[:,0]+risetime, tail[:,1], tail[:,2], fmt='g.')
#plt.plot(t, ft.edp_nomc(t, t0), label='t0='+str(t0))
#plt.plot(tail[:,0]+risetime, ft.edp_nomc(tail[:,0]+risetime, 28), label='t0=28')



plt.legend(loc=0, numpoints=1)
plt.xlabel('Time since explosion (days)')
plt.ylabel('Bolometric Flux ($\cdot 10^{43}$ erg s$^{-1}$)')


#chisq plot
#plt.figure(2); tarr=np.linspace(0, 50, 100)
#plt.plot(tarr, np.array([chisq_by_hand(ft.edp_nomc,i, tail, err=True) for i in tarr]))

#chisq plot (optional, set this argument to True)

plotchi = False
if plotchi:
    plt.figure(2)
    tarr=np.linspace(0, 50, 100)
    plt.plot(tarr, np.array([chisq_by_hand(ft.edp_nomc,i, tail) for i in tarr]))
plt.ylim(0, 1.5)
plt.xlim(0, 150)
plt.savefig("tailfit_transparency_08hv.pdf")
plt.show()
