"""
!!!! Requires pack for line fitting !!!!
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from pack import bol

def plt_arr_crt(arr, num):
	r = np.array([float(arr[arr[:,0]==i[0]][0][num]) for i in mej if i[0] in mej])
	return r

mej1 = np.loadtxt('../ejecmass_bolerr.txt', dtype='string', usecols=(0, -2, -1))
mej = mej1[1:] ##improper errors on SN2005eq, hence, remove

sbv = np.loadtxt('../sbv_all_b14.txt', dtype='string')

sbv_arr = plt_arr_crt(sbv, 1)
esbv_arr = plt_arr_crt(sbv, 2)
mej_arr = plt_arr_crt(mej, 1)
emej_arr = plt_arr_crt(mej, 2)


vs = np.vstack([sbv_arr, esbv_arr, mej_arr, emej_arr]).T
rn = np.random.normal

vals, errs = bol.linfits().freq_odr(vs)

def calc_mej_fromsbv(sbv):
	a, b = sbv
	ar = np.array([rn(vals[0], errs[0])*rn(a, b) + rn(vals[1], errs[1]) for k in range(100000)])	
	return np.mean(ar), np.std(ar)

sne_sublum = np.loadtxt('all91bg.txt', dtype='string')
sbv_sublum = np.array([[float(i[1]), float(i[2])] for i in sbv if i[0] in sne_sublum])

mej_sublum = np.array([calc_mej_fromsbv(i) for i in sbv_sublum])
mej_tot = np.array([calc_mej_fromsbv([float(i[1]), float(i[2])]) for i in sbv])

print mej_sublum
#plt.hist(mej[:,1].astype('float32'), alpha=0.5, bins=np.arange(0.4, 1.8, 0.2))
plt.hist(mej_sublum[:,0], alpha=0.5, bins=np.arange(0.4, 1.8, 0.1))
plt.hist(mej_tot[:,0], alpha=0.5, bins=np.arange(0.4, 1.8, 0.1))
plt.show() 







