#Variables for main.py
from constants import m_H

#Star-specific variables
stellar_tracks = 'Hidalgo' #Stellar luminosity reference options are 'Baraffe' (Baraffe et al., 2015) and 'Hidalgo' (Hidalgo et al., 2018)
norm_lum = True   #Normalize interpolated luminosity to match stellar luminosity
m_star = 0.437    #solar masses
r_star = 0.4232   #solar radii
logg = 4.826      #surface gravity, cm/s2
Teff = 3573       #K
FeH = -0.08       #stellar metallicity [dex]
J_mag_star = 9.706#2MASS J magnitude
d_star = 36.011   #parsecs

#Planet and atmosphere-specific variables
a_planet = 0.1037 #Semi-major axis, au
albedo = 0.0       #Bond albedo (dimensionless)

#planet core and envelope mass fractions, envelope composition
#a present planet w/ 20% envelope @ 90% H2O, 10% H2 evolves from a planet with an envelope fraction of 2.62447e-01 and 
#     composition = [0.31764360070863473, 4.386069211982582e-51, 0.6823563992913653, 2.958440455988615e-51]; initial R = 1.698912 and M = 2.06086
# 99.5% core_frac, all H2 - initial envelope fraction = 4.50385e-02 and composition = [1.0, 0.0, 0.0, 0.0]; initial R = 1.688117 and M = 2.019254
core_frac = 0.9   #core mass fraction - core and envelope fractions are now computed on the fly from bulk composition
core_den = 5.E3   #core mean density, kg/m3
envelope_frac = 1-core_frac #do the envelope by subtraction - core and envelope fractions are now computed on the fly from bulk composition
envelope_species= ['H2', 'He',  'H2O', 'O', 'CO2'] #component molecular masses [amu]
envelope_compm= [2, 4, 18, 16, 44] #component molecular masses [amu]
envelope_compH= [1, 0., 2/3, 0., 0.] #fraction of atoms that are H in each species

#Variables related to escape
do_emp_sat = True   #empirically match the saturated and subsaturated regime rates from Sanz-Forcada et al. (2011), eq. 5
do_euv_sat = True   #Use the empirical EUV flux saturation from Peacock et al. (2020) (HAZMAT VI)
do_emp_scale = True  #Use the empirical X and EUV flux evolution from Peacock et al. (2020) (HAZMAT VI) instead of Sanz-Forcada et al.
xuv_threshold = 1.E4 #erg/cm2/s = flux threshold for X-ray driven escape fluxes; 
#                       notional value ~7.E3 from Owen & Jackson (2012; MNRAS)
p_photo = 2.E3   #20 mbar = 2E-2 bar * 1E5 = 2E3 Pa transit 'radius' for planet
p_xuvcofac = m_H/1.89E-14 #kg/m2; non-varying cofactor for P_base (XUV photosphere pressure) from Lopez et al (2017; MNRAS)
variable_efficiency = False #consider if the mass loss efficiencies vary with other parameters
hnu = 1.025E-18 #J, equivalent of 6.4 eV (~195 nm) (Owen & Alvarez, 2016); mean XUV photon energy that is responsible for heating the atmosphere

#numerical variables
mode = 'forward'#'reverse'  #Use the escape fluxes in reverse mode to 'add' mass to present-day planet, or 'forward' to evolve some initial state
