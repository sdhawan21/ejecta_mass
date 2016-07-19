#plot observed and model Mej-Mni
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

from glob import glob
home = os.path.expanduser("~")
pt = home+'/workspaces/source_nimass/fast_decl/nimass_ejmass/redfree_estimate/'
fig = plt.figure(figsize=(5.5,6))

models = glob('models/*.dat')
models = np.delete(models, 3)
models = np.delete(models, 3)
labels = ['Blondin+2013:DDC', 'Fink+2014:DEFL', 'Woosley+2011:DoubleDet', 'Fink+2010:DoubleDet','Sim+2010:sub-mc Det']
data = np.loadtxt(pt+'ejecta_realerr_fin.txt', usecols=(2,3,4,5))
forms = ['rs', 'gD', 'yo', 'bD', 'mx']
print models
for i, mod in enumerate(models):
	f = forms[i]
	mval = np.loadtxt(mod, skiprows=1)
	plt.plot(mval[:,0], mval[:,1], f, label=labels[i])
plt.errorbar(data[:,-2], data[:,0], xerr=data[:,-1], yerr=data[:,1], fmt='k^', label='Calculated from data')
plt.xlabel("$\mathrm{M_{^{56}Ni}}$")
plt.ylabel("$\mathrm{M_{ej}}$")
plt.legend(loc=0, numpoints=1, prop={'size':10})
plt.savefig("ejectamass_nickel_modelcomp.pdf")
