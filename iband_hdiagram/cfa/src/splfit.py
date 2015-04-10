import numpy as np
import easygui as eg

from scipy.interpolate import interp1d
from glob import glob
def main():
	indir=glob('../data/*.ascii')
	z=np.loadtxt('../cfa3_redshift.ascii', dtype='string')
	imaxarr=[]
	for i in indir:
		ph, mag, magerr= np.loadtxt(i, unpack=True)
		pp=interp1d(ph, mag)
		l=np.linspace(ph.min(), ph.max())
		inmag=pp(l)
		sn=i[8:-8]
		
		redsn=float(z[z[:,0]==sn][0][1])
		if min(inmag) != inmag[0] and redsn > 0.01 and min(ph) < -3:
			imaxarr.append([sn, redsn, min(inmag)])
	np.savetxt('../imax_samp1.ascii', imaxarr, fmt="%s")
	print imaxarr, len(imaxarr)
main()