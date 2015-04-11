"""
Source code to fit splines to first max

Dependencies: easygui
"""
import easygui as eg
import matplotlib.pyplot as plt
import numpy as np


from scipy.interpolate import interp1d
from glob import glob
def main():
	indir=glob('../data/*.ascii')
	z=np.loadtxt('../cfa3_redshift.ascii', dtype='string')
	imaxarr=[]
	for i in range(len(indir)):#indir:
		ph, mag, magerr= np.loadtxt(indir[i], unpack=True)
		
		sn=indir[i][8:-8]
		"""
                #make selection using easygui (need to learn better gui package)
                ranarr=eg.multchoicebox('choose a value for the fit',sn, ph)
                ranarr=np.array(ranarr)
                sortarr=sorted(ranarr.astype('float32'))
               
                #use gui selection
                ph1=ph[(ph > sortarr[0] ) & (ph < sortarr[1])]
                mag1=mag[(ph > sortarr[0] ) & (ph < sortarr[1])]
                """
                #interpolate with scipy
                pp=interp1d(ph, mag)
		l=np.linspace(ph.min(), ph.max(), 1000)
		inmag=pp(l)

                plt.plot(ph, mag, 'go')
                plt.plot(l, inmag, 'k:')
		plt.ylim(plt.ylim[::-1])
                #plt.show()
		redsn=float(z[z[:,0]==sn][0][1])
		if min(inmag) != inmag[0] and redsn > 0.01:
			imaxarr.append([sn, redsn, min(inmag)])
		plt.savefig('../img/'+sn+'.png')
	np.savetxt('../imax_samp1.ascii', imaxarr, fmt="%s")
	print imaxarr, len(imaxarr)
main()
