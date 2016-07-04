import numpy as np
import matplotlib.pyplot as plt
import sys
import os

from scipy.optimize import curve_fit
from pack import bol

bp = bol.bol_func().bolpeak
class energy_dep:
	
	def __init__(self):
		self.eni = 3.70e10	#erg/g
		self.eco = 6.76e9	#erg/g
		self.lni = 1.3152e-6
		self.lco = 1.0325e-7

	def radius(self, t, vexp):
		return vexp*t

	def edep_pos(self, t, mni):
		a = self.eco*(np.exp(-self.lco*t)-np.exp(-self.lni*t))
		return 0.08*mni*a

	def edep_gam(self, t, mni):
		sni=self.eni*np.exp(-self.lni*t)
		sco = self.eco*(np.exp(-self.lco*t)-np.exp(-self.lni*t))
		return mni*(sni+0.82*sco)
	
	def f(self, n, x0=0.15):
		a = x0**3/3
		b = pow(x0, n)*(1 - pow(x0, 3-n))/(3-n)
		return a+b
		
	def rho(self, t, mej, n, vexp, x0=0.15):
		coef=mej/(4*np.pi*self.f(n))	
		return coef/pow(self.radius(t, vexp), 3)
		
	def tau_pos(self, t, kpos, mej, n, vexp, x0=0.15):
		e = 1+((1 - pow(x0,n-1))/(n-1))
		return e*self.rho(t,mej,n,vexp)*kpos*self.radius(t, vexp)*x0
	
	def tau_gam(self, t, kgam, mej, n, vexp, x0=0.15):
		e = 1+((1 - pow(x0,n-1))/(n-1))
		return e*self.rho(t,mej,n,vexp)*kgam*self.radius(t, vexp)*x0
	
	def edep_tot(self, t, kgam,mej,n,mni,kpos,vexp=30000):
		g = self.edep_gam(t,mni)*(1-np.exp(-self.tau_gam(t, kgam, mej,n,vexp)))
		
		p = self.edep_pos(t, mni)*(1-np.exp(-self.tau_pos(t, kpos, mej,n,vexp)))
		return g+p

def main():	
	fname = sys.argv[1]	
	rt = float(sys.argv[2])
	
	bollc = np.loadtxt(fname)
	mmax, tmax = bp(fname)
	bollc[:,0]-=tmax
	bollc[:,0]+=rt
	
	tail = bollc[bollc[:,0]>40.]
	popt, pcov = curve_fit(energy_dep().edep_tot, tail[:,0], tail[:,1], sigma=tail[:,2], p0=[0.03, 1., 2., .5, 1.0])
	print popt
	
	t=np.arange(bollc[:,0].min(), bollc[:,0].max(), 1)
	edep = energy_dep().edep_tot(t, 0.03, 1.0,2.1,0.6,1.0)
	
	#plt.plot(bollc[:,0], bollc[:,1], 'bs')
	plt.plot(t, edep)
	plt.show()
if __name__=="__main__":
	main()

