#Variables for main.py
from constants import m_H

#Star-specific variables
stellar_tracks = 'Hidalgo' #Stellar luminosity reference options are 'Baraffe' (Baraffe et al., 2015) and 'Hidalgo' (Hidalgo et al., 2018)
norm_lum = True   #Normalize interpolated luminosity to match stellar luminosity
m_star = 0.437    #solar masses
r_star = 0.4232   #solar radii
age_star = 7.9    #stellar age in Gyr
l_star = 0.02629  #solar luminosities
logg = 4.826      #surface gravity, cm/s2
Teff = 3573       #K
FeH = -0.08       #stellar metallicity [dex]
J_mag_star = 9.706#2MASS J magnitude
d_star = 36.011   #parsecs

#Planet and atmosphere-specific variables
m_planet = 1.9#1.998816 #1.9    #Earth masses
r_planet = 1.67#1.682758 #1.67   #Earth radii
a_planet = 0.1037 #Semi-major axis, au
albedo = 0.0       #Bond albedo (dimensionless)

#planet core and envelope mass fractions, envelope composition
#a present planet w/ 10% envelope @ 90% H2O, 10% H2 evolves from a planet with 
#     an envelope of 14.4493% and composition 0.31347321438588577 H2 and 0.6865267856141143 H2O, assuming R = 1.682758 and M = 1.998816
core_frac = 0.9#1-0.14493   #core mass fraction
core_den = 5.E3   #core mean density, kg/m3
envelope_frac = 1-core_frac #do the envelope by subtraction
envelope_comp = [0.1,0,0.9,0] #[0.31347321438588577, 0.0, 0.6865267856141143, 0.]  #surface mixing ratios - must sum to 1
envelope_species= ['H2', 'He',  'H2O', 'O'] #component molecular masses [amu]
envelope_compm= [2, 4, 18, 16] #component molecular masses [amu]
envelope_compH= [1, 0., 2/3, 0.] #fraction of atoms that are H in each species
efficiencies = [0.1, 0.1, 0.01, 0.01] #escape efficiency epsilon for each component

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
mode = 'reverse'  #Use the escape fluxes in reverse mode to 'add' mass to present-day planet, or 'forward' to evolve some initial state
