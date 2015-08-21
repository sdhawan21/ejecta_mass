import numpy as np
import matplotlib.pyplot as plt
import shutil
import sys

from scipy.stats import pearsonr, ks_2samp
from pack import bol

def ir_frac(sn):
	#fraction of flux emitted in the NIR (YJH) as a function of phase
	bp = bol.bol_func().bolpeak
	tot='_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
	pt = '/home/sdhawan/bol_ni_ej/lcbol_distrib/'
	
	u = pt+sn+tot
	
	ir = 	tot='_lcbol_Y_CSPJ_CSPH_CSP.dat'
	i = pt+sn+ir
	
	lc = np.loadtxt(u) - bp(u)[1]
	
	irlc = np.loadtxt(i)-bp(u)[1]
	fracarr=[]
	pharr =[]
	for i in lc:
		if i[0] in irlc[:,0]:
			fracarr.append(irlc[irlc[:,0]==i[0]][0][1]/i[1])
			pharr.append(i[0])
	return 	np.array(pharr), np.array(fracarr)

def ir_at_max(sn):
	p, f = ir_frac(sn)
	f1=f[abs(p) == min(abs(p))][0]
	return f1
def main():
	pt = '/home/sdhawan/bol_ni_ej/'
	
	sn =sys.argv[1]
	p, f=ir_frac(sn)
	print "the NIR fraction at bolometric maximum is:", round(ir_at_max(sn)*100, 2), "%"

	plt.plot(p, f)
	plt.show()
	return 0
	ir = np.loadtxt('../ejecmass.txt', usecols=(-2, -1))
	noir = np.loadtxt('../ejecmass_noir.txt', usecols=(-2, -1))
	
	print ks_2samp(ir[:,0], noir[:,0])
	plt.hist(ir[:,0], histtype='step')
	plt.hist(noir[:,0], histtype='step')
	#plt.show()
	
	
if __name__ == "__main__":
	main()
