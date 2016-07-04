"""
Use the Bateman equation to evaluate the contribution of less abundance radioactive nuclei to the late time light curve.
<<<<<<< HEAD
=======
--> Seitenzahl+ 2009 has the decay energies for the isotopes
--> Seitenzahl+ 2014 has the e-folding timescales

--> Graur+2016 presents the equation 

>>>>>>> 16bd5f60542a4bb14357e4a409caafb34c2a1020
"""
import numpy as np
import matplotlib.pyplot as plt

co57_kev = 7.594 + 10.22 + 3.598
co56_kev = 3.355 + 0.374 + 115.7 + 1.588
fe55_kev = 3.973 + 1.635


<<<<<<< HEAD

a_co = 27
a_fe = 26

m_co57 = 0.05
m_co56 = 0.6
m_fe55 = 0.03
=======
#atomic numbers
a_co = 27
a_fe = 26

#'test' values for abundances
m_co57 = 0.02
m_co56 = 0.8
m_fe55 = 0.01
>>>>>>> 16bd5f60542a4bb14357e4a409caafb34c2a1020

lam_co7 = 2.551e-3
lam_co6 = 8.975e-3
lam_fe  = 6.916e-4

<<<<<<< HEAD
=======
#from Graur+2016
>>>>>>> 16bd5f60542a4bb14357e4a409caafb34c2a1020
def l_bateman(t, mass, dec_energy, lambda_a, a, B=0.235):
	 fac = 2.221*B/a
	 return fac*lambda_a*dec_energy*np.exp(-lambda_a*t)*1e43
	 
t_axis = np.linspace(600, 1200)
l_co57 = l_bateman(t_axis, m_co57, co57_kev, lam_co7, a_co)
l_co56 = l_bateman(t_axis, m_co56, co56_kev, lam_co6, a_co)
l_fe55 = l_bateman(t_axis, m_fe55, fe55_kev, lam_fe, a_fe)

tot = l_co57 + l_co56 + l_fe55

<<<<<<< HEAD
plt.plot(t_axis, l_co57)
=======
plt.plot(t_axis, l_co57, label="57Co "+str(m_co57))
plt.plot(t_axis, l_co56, label="56Co "+str(m_co56))
plt.plot(t_axis, l_fe55, label="55Fe "+str(m_fe55))
plt.plot(t_axis, tot, label="All isotopes")

plt.legend(loc=0, numpoints=1)
>>>>>>> 16bd5f60542a4bb14357e4a409caafb34c2a1020
plt.show()
