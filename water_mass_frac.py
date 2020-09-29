#This code will scan over parameter space to retrieve values for the core and water
# fractions based on Noack et al., 2016. 

from math import pi,log10,exp
from constants import G,kb,m_Earth,r_Earth,m_H
from modules import bisect2

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
            if i > 1000:
                  #guess = -1000
                  break
      return guess,function(guess)

# m_planet = [0,0.6,1.9,4.2,6.4] #M_Earth -2, -1, 0, 1, 2 sigma
# r_planet = [1.563,1.673,1.76] #R_Earth 2 sigma
m_planet = 1.9#6.4 
r_planet = 1.673#1.563
Rp = r_planet*r_Earth
Mp = m_planet*m_Earth
g0 = G*Mp/(Rp**2.)
bulk_den = (Mp/(4*pi*Rp**3./3))/1E3 #convert to g/cm3
print(f'Planet radius = {Rp:10.2e} ({r_planet:4.2f} Earth radii); mass = {Mp:10.2e} ({m_planet:4.2f} Earth masses); surface gravity = {g0:8.2f} m/s2.')
print(f'Planet bulk density = {bulk_den:5.2f} g/cm3 (for reference, Earth is ~5.5 g/cm3).')
mu = 2. #Hydrogen envelope 
T_eq = 500. #TOI-1266c's initial equlibrium temperature, before star enters MS
R_B = (G*Mp/(2*kb*T_eq/(m_H*mu)))/r_Earth #R_B = G*M/2*cs**2
print(f"Assuming that mu = {mu:3.1f}, for a T_eq = {T_eq:6.2f} K, the Bondi radius = {R_B:6.3f} R_earth")

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

##Zeng et al., 2013; 2016 - Fe and Earth-like compositions
##CMF_Z = (1/0.21)*(1.07 - r_planet/m_planet**(1/3.7)) #<-- this is a negative number for TOI-1266c's reported M/R

Rfe_MG = lambda mp: 0.815*(mp)**(1/4.176) #100% iron
Rrock_MG = lambda mp: 1.007*(mp)**(1/3.7) #100% silicate
Rrock1_MG = lambda mp: 1.41*(mp)**(1/3.905) #1% H2 on rocky core
RH2O_MG = lambda mp: 1.41*(mp)**(1/3.905) #100% H2O
RH2_MG = lambda mp: 4.106*(mp)**(1/5.01) #100% H2 planet 
Rremnant_MG = lambda mp: 0.469*(mp)**(1/3) #From Mocquet et al., 2014, for super-dense remnants
#Modirrousta-Galian et al., 2020 - The Bimodal Distribution of exoplanet radii - has parameterizations from Zeng et al. (2013; 2016)
print(f"Modirrousta-Galian radii fits to Zeng et al. curves:"
      + f"\n       100% Fe       = {Rfe_MG(m_planet):6.4f} R_earth"
      + f"\n       100% silicate = {Rrock_MG(m_planet):6.4f} R_earth"
      + f"\n           + 1% H2   = {Rrock1_MG(m_planet):6.4f} R_earth"
      + f"\n       100% H2O      = {RH2O_MG(m_planet):6.4f} R_earth"
      + f"\n       100% H2       = {RH2_MG(m_planet):6.4f} R_earth"
      + f"\n       Remnant core  = {Rremnant_MG(m_planet):6.4f} R_earth")

AN = lambda xfe,xh2o: 7121 - 20.21*xfe + xh2o*(15.23 + 0.239*xfe)
ANFE = lambda xfe: 7121 - 20.21*xfe + (1-xfe)*(15.23 + 0.239*xfe)
ANH2O = lambda xh2o: 7121 + xh2o*(15.23)
CN = lambda xh2o: 0.2645 + 0.00048*xh2o
f_RpN = lambda xfe,xh2o,mass: 1000*AN(xfe,xh2o)*mass**CN(xh2o)/r_Earth
#Noack et al. (2016) - Water-rich planets: How habitable... Eq. 4 and Table 4
#     Rp normally in km, Mp in Earth masses
fracFe = 55
fracSi = 45
fracH2O = 100-fracSi-fracFe
print(f'Noack et al. (2016) prediction for {fracFe:2.0f}% Fe, {fracSi:2.0f}% silicate, and {fracH2O:2.0f}% H2O = {f_RpN(fracFe,fracH2O,m_planet):5.2f} Earth radii')
bFeN,bH2ON,brpN = bisect2(0,100,0,100,'Fe',r_planet,m_planet,f_RpN,1.E-4)
print(f'    From bisect: {bFeN:2.0f}% Fe, {100-bFeN-bH2ON:2.0f}% silicate, {bH2ON:2.0f}% H2O -> {brpN:5.2f} Earth radii')

exit()
#Duan & Zhang, 2006 - Equation of state of the H2O, CO2, and H2O–CO2 systems up to 10 GPa and 2573.15 K
#Added alpha, beta, and gamma as last three terms in 'a'
a = [4.68071541E-2, -2.81275941E-1, -2.43926365E-1, 1.10016958E-2, -3.86603525E-2, 9.30095461E-2, -1.15747171E-5,
      4.19873848E-4, -5.82739501E-4, 1.00936000E-6, -1.01713593E-5, 1.63934213E-5, -4.49505919E-2, -3.15028174E-1, 1.25E-2]
Rweird = 83.14467 #cm3 bar/K /mol
Vc = Rweird*647.25/221.19 #bar?

T_surf = 2E3
V = 1E6 # cm3 equivalent of 1 m3
P = 1.4E4 #bar equivalent to 1.4E9 Pa

Tr = T_surf / 647.25 #K/K
B = lambda T: a[0] + a[1]/T**2 + a[2]/T**3
C = lambda T: a[3] + a[4]/T**2 + a[5]/T**3
D = lambda T: a[6] + a[7]/T**2 + a[8]/T**3
E = lambda T: a[9] + a[10]/T**2 + a[11]/T**3
F = lambda T: a[12]/T**3
#Z = PV/RT
Z = lambda V: (1 + B(Tr)*Vc/V + C(Tr)*Vc**2/V**2 + D(Tr)*Vc**4/V**4 + E(Tr)*Vc**5/V**5 
                  + (F(Tr)*Vc**2)*(a[13] + a[14]*Vc**2/V**2)*exp(-1*a[14]*Vc**2/V**2) )
# Z(V)*Rweird*T_surf/P == V (must equal the input V!
#between 52 and 53 = 52.4112 cm3/mol for P=1.4E4 bar at 2,000 K

print(f"Pressure test = {[Z((v+5241120)/100000)*Rweird*T_surf/P for v in range(0,10)]}")
