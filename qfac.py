<<<<<<< HEAD

import numpy as np
from scipy.integrate import quad, tplquad

velocity_range = [0., np.inf]
=======
import numpy as np

from scipy.integrate import quad, tplquad

#velocity range in units of e-folding velocity 
ve = 2700.
v_outer = 15000.
velocity_range = [0.0, np.inf] #v_outer/ve]

>>>>>>> 16bd5f60542a4bb14357e4a409caafb34c2a1020
a = velocity_range[0]; b=velocity_range[1]

q_formfac_num = lambda zs,u,z: z**2*np.exp(-z)*np.exp(-np.sqrt(z**2+zs**2+2*u*zs*z))
	
zs1,zs2 = a , b
u1, u2 = lambda zs:-1, lambda zs:1
z1, z2 = lambda zs,u: a,  lambda zs,u: b

tplint = tplquad(q_formfac_num, zs1, zs2, u1, u2, z1, z2)
print tplint

q_formfac_den = lambda z: z**2*np.exp(-z)
den_int = quad(q_formfac_den, a, b)

qval = tplint[0]/(den_int[0]*2)
print qval
