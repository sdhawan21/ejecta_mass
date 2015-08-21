#calculate ejecta masses for model bolometric light curves
import numpy as np
import matplotlib.pyplot as plt
import sys


from pack import bol, fid_time


def main():
	mod = sys.argv[1]
        fwrite = sys.argv[2]
        ir = sys.argv[4]
        
       	if ir == '1':
       		suf = '_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'
       	else:
       		suf = '_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSP.dat'
       
	f = '/home/sdhawan/bol_ni_ej/lcbol_distrib/'+mod+suf

	bp = bol.bol_func().bolpeak
	m,t=bp(f)
        
        
        dfile = np.loadtxt('/home/sdhawan/tests_paper/ni/files_snpy/tmax_dm15.dat', dtype='string')
        
        fout= open('model_recons.txt', 'a')

	rt = t
	mni = m/(2e43)
       
	ft = fid_time.fid_time(mni)
	
	#use a signature for observed data
	if 'SN' in mod:
		dm = float(dfile[dfile[:,0] == mod][0][3])
		rt = 16.5- 5*(dm -1.1)
	
	print "The Nickel mass is:", mni
        print "The rise time is :", rt
	
	bb = ft.val_calc(mod, '/home/sdhawan/', rt, ir)
        
        if fwrite == True:
                fout.write(mod+'\t'+str(bb[1])+'\n')
        fout.close()
        print suf
        print bb, t
        
        plotlc = sys.argv[3]
       
	
	if plotlc == True:
	
		sn = np.loadtxt(f)
		ph = sn[:,0] - bp(f)[1]
		if len(sn[0]) == 2:
			plt.plot(ph, sn[:,1])
		elif len(sn[0]) == 3:
			plt.errorbar(ph, sn[:,1],  sn[:,2], fmt='go')
		else:
			print "What is this file size, recheck: Stopping the script"
			return 0
		plt.plot(ph, ft.edp_nomc(sn[:,0]-t, bb[0]))
		
			
	
		plt.yscale('log')
		plt.show()
	
	
main()




