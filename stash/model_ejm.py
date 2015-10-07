import numpy as np
import sys


from pack import bol, fid_time


def main():
	mod = sys.argv[1]
        fwrite = sys.argv[2]
	f = '/Users/lapguest/bol_ni_ej/lcbol_distrib/'+mod+'_lcbol_u_CSPB_CSPV_CSPr_CSPi_CSPJ_CSPH_CSP.dat'

	bp = bol.bol_func().bolpeak
	m,t=bp(f)
        
        
        fout= open('model_recons.txt', 'a')

	rt = t
	mni = m/(2e43)
        print "The Nickel mass is:", mni
        print "The rise time is :", rt
	ft = fid_time.fid_time(mni)

	bb = ft.val_calc(mod, '/Users/lapguest/', rt)
        
        if fwrite:
                fout.write(mod+'\t'+str(bb[1])+'\n')
        fout.close()
        print bb, t
        

main()




