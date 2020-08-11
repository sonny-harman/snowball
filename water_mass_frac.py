#This code will scan over parameter space to retrieve values for the core and water
# fractions based on Noack et al., 2016. 

from math import pi,log10
from planet import r_planet,m_planet
from constants import G,m_Earth,r_Earth,m_H

def bisect(low,high,match,function,threshold):
      guess = (low + high)/2
      eps = abs(match - function(guess))
      i = 0
      while eps > threshold:
            i += 1
            if function(guess) > match:
                  high = guess
            if function(guess) < match:
                  low = guess
            guess = (low + high)/2
            eps = abs(match - function(guess))
            if i > 100:
                  guess = -1000
                  break
      return guess,function(guess)

Rp = r_planet*r_Earth
Mp = m_planet*m_Earth
g0 = G*Mp/(Rp**2.)
bulk_den = (Mp/(4*pi*Rp**3./3))/1E3 #convert to g/cm3
print(f'Planet radius = {Rp:10.2e} m ({m_planet:4.2f} Earth radii); surface gravity = {g0:8.2f} m/s2.')
print(f'Planet bulk density = {bulk_den:5.2f} g/cm3 (for reference, Earth is ~5.5 g/cm3).')

#Fortney et al. (2007) - Planet radii across five orders of magnitude - Eq. 7 (note the Erratum correction!)
ice_low = 0.; ice_high = 1.; ice_mf = (ice_low + ice_high)/2
lm = log10(m_planet)
Rf = lambda i: (0.0912*i + 0.1603)*lm**2. + (0.3330*i + 0.7387)*lm + (0.4639*i + 1.1193)
Rf_fit = Rf(ice_mf)
ice_mf,Rf_fit = bisect(ice_low,ice_high,r_planet,Rf,1.E-3) 
print(f"Fortney et al. (2007) radius prediction for {ice_mf*100:4.2f}% water ice = {Rf_fit:5.2f} Earth radii")

#Fu et al. (2009) - The interior dynamics of water planets - Eq. 17 and Table 2
#     Doesn't reach the observed radius, but important for comparison
water_mf = [0.02, 25, 50] #Percent water by mass
Af = [0.994, 1.157, 1.255]
Bf = [0.266, 0.253, 0.251]
Rpf = [a*m_planet**b for a,b in zip(Af,Bf)]
print("Fu et al. (2009) radius prediction for varying water mass fractions:")
for i in range(len(Rpf)):
      print(f"\t {water_mf[i]:4.2f}% H2O = {Rpf[i]:5.2f} Earth radii")


#Nettelman et al., 2010 - Interior models of GJ436b - Eq. 1 and coeff. in text
#     (0.25 M_E rocky core, T set at 1 bar)
N436_iso = [1.6586, 0.9950, 0.1549, 0.]         #T = 1000 K @ 1 bar, isothermal
N436_ad = [2.8210, -0.2928, 0.9037, -0.176]     #T = 1000 K @ 1 bar, fully adiabatic
Rni = sum([N436_iso[i]*log10(m_planet)**i for i in range(len(N436_iso))])
print(f"Nettelman et al. (2010) radius prediction for 1000K isothermal atm. = {Rni:5.2f} Earth radii")
Rna = sum([N436_ad[i]*log10(m_planet)**i for i in range(len(N436_ad))])
print(f"Nettelman et al. (2010) radius prediction for 1000K adiabatic atm. = {Rna:5.2f} Earth radii")

#Noack et al. (2016) - Water-rich planets: How habitable... Eq. 4 and Table 4
#     Rp normally in km, Mp in Earth masses
XFe = 0.; XH2O = 77.
AN = lambda xfe,xh2o: 7121 - 20.21*xfe + xh2o*(15.23 + 0.239*xfe)
ANFE = lambda xfe: 7121 - 20.21*xfe + (1-xfe)*(15.23 + 0.239*xfe)
ANH2O = lambda xh2o: 7121 + xh2o*(15.23)
CN = lambda xh2o: 0.2645 + 0.00048*xh2o
f_RpN = lambda xfe,xh2o: 1000*AN(xfe,xh2o)*m_planet**CN(xh2o)/r_Earth
print(f'Noack et al. (2016) prediction for 0% Fe and 100% H2O = {f_RpN(0,100):5.2f} Earth radii')
