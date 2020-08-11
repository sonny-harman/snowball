from math import pi

#Physical constants 
m_Earth = 5.9724E24     #kg
r_Earth = 6.378E6       #m
m_Sun = 1.9891E30       #kg
r_Sun = 6.963E8         #m
l_Sun = 3.839E33        #!!! erg/s !!!
G = 6.674E-11           #m3/kg/s2
kb = 1.38064E-23        #kg*m2/s2/K
Na = 6.022E23           #molec/mol
Rconst = kb*Na          #J/K/mol
h = 6.626E-34           #J*s (kg*m2/s2)
c = 2.99E8              #m/s
sigma = 5.67E-8         #W/m2/K4

m_H = 1.67E-27          #kg
#ev = 1.60218e-19      #1 ev = 1.60218e-19 Joule

#conversion factors
au2m = 1.495978707E+11  #convert au to m
m2cm = 100.             #convert from m to cm
cm2m = 1/m2cm           #convert from cm to m
yr2s = 365*24*3600      #convert year to seconds
s2yr = 1/yr2s           #convert seconds to year
ergcm2s2Wm2 = 1E-3      #convert from ergs/cm2/s to W/m2

#constant vales
flux_Earth = l_Sun/(4*pi*(au2m*m2cm)**2.)

#also include uncertainties from Sanz-Forcada et al. (2011)
euv_a = 0. #+/-0.11; first term in eq. 5
euv_b = 0. #+/-0.15; second term in eq. 5
