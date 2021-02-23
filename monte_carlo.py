from math import log10
from numpy import logspace,linspace
from main import run_escape
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from smt.sampling_methods import LHS
import constants

#TODO replace the following variables in planet.py with a reference to the ones generated here
#? stellar_tracks = ['Hidalgo', 'Baraffe']
# age_star = [2.7,7.9,12.1] #Gyr 2 sigma
# l_star = [0.02554,0.02629,0.0270] #L_sun 2 sigma
# m_planet = [0,0.6,1.9,4.2,6.4] #M_Earth -2, -1, 0, 1, 2 sigma
# r_planet = [1.563,1.673,1.76] #R_Earth 2 sigma
# core_frac, envelope_frac - set by max. H2O given R,M
# envelope_comp - H2,He,H2O,O,CO2 - no O to start, vary H2O, H2, CO2, set He as x*H2?
# efficiencies = [0.01-0.3, 0.1, 0.01-0.3, 0.01, 0.01] - only test ranges for H2O, H2

#TODO save variables from main.py (turned into a function call at the end of main?)
#initial and final states for mass, radius, envelope composition; total H lost
#stellar age, luminosity

running_MC = True
ntimes = 10000

#Age, Luminosity, Mass, Radius, H2, He, H2O, CO2, H2 escape eff., H2O escape eff., other eff.
limits = np.array([[2.7,12.1],[0.02554,0.0270],[1.575,6.4],[1.563,1.76],[-6,0],[-6,0],[-6,0],[-6,0],[0.01,0.4],[0.01,0.4],[0.01,0.1]])
sampling = LHS(xlimits=limits)
rand_vals = sampling(ntimes)

age = []; lum=[]; mass = []; radius = []; envelope = []; eff = []
mass_i = []; mass_c = []; mass_f = []
radius_i = []; radius_c = []; radius_f = []
core_i = []; core_c = []; core_f = [];
env_i = []; env_c = []; env_f = [];
env_comp_i = []; env_comp_c = []; env_comp_f = [];
dmass = []

for i in range(ntimes):
      #Transform gas mass fractions from log to linear
      for k in range(4,8):
            rand_vals[i,k] = 10**(rand_vals[i,k])
      #Solar helium to hydrogen is 7.8369e-02/9.2076e-01 = 8.51%, from Asplund (2005); 73.8% by mass H, 24.85% by mass He
      He_abundance = 0.3367*rand_vals[i,4]*rand_vals[i,5] #ind=4 is the hydrogen mass fraction, so this gives us a He abundance anywhere
      #     from 0% to solar He/H (when ind=5 == 1).
      #Now, renormalize the gas composition
      tot_gases = rand_vals[i,4] + He_abundance + rand_vals[i,6] + rand_vals[i,7] 
      He_abundance = He_abundance/tot_gases
      rand_vals[i,4] = rand_vals[i,4]/tot_gases
      rand_vals[i,5] = He_abundance
      rand_vals[i,6] = rand_vals[i,6]/tot_gases
      rand_vals[i,7] = rand_vals[i,7]/tot_gases
      comp = rand_vals[i,4:8]
      comp = np.insert(comp,3,0)
      print(f"Run {i:4} for m_planet = {rand_vals[i,2]:6.2f} and atmospheric comp. = {comp}")
      #print(rand_vals[i,:])
      other_effs = rand_vals[i,10]
      initial_comp = [rand_vals[i,4],rand_vals[i,5],rand_vals[i,6],0,rand_vals[i,7]]
      total_effs = [rand_vals[i,8],other_effs,rand_vals[i,9],other_effs,other_effs]

      mi,mc,mf,ri,rc,rf,ci,cc,cf,ei,ec,ef,eci,ecc,ecf,mass_change = run_escape(
                        rand_vals[i,0],rand_vals[i,1],rand_vals[i,2],rand_vals[i,3],
                        initial_comp,total_effs)
                        #MC_age[i],MC_lum[i],MC_mass[i],MC_radius[i],MC_envelope[i],MC_eff[i])
      
      if mi != rand_vals[i,2]:
            print('That is not good.',mi,rand_vals[i,2])

      age.append(rand_vals[i,0]); lum.append(rand_vals[i,1])
      eff.append(total_effs)
      mass_i.append(mi); mass_c.append(mc); mass_f.append(mf)
      radius_i.append(ri); radius_c.append(rc); radius_f.append(rf)
      core_i.append(ci); core_c.append(cc); core_f.append(cf)
      env_i.append(ei); env_c.append(ec); env_f.append(ef)
      env_comp_i.append(eci); env_comp_c.append(ecc); env_comp_f.append(ecf)
      dmass.append(mass_change)

df = pd.DataFrame({'Current Stellar Age [Gyr]':age,
      'Luminosity [L/L$_{\odot}$]':lum,
      'Loss Efficiency':eff,
      'Initial Mass [M$_{\oplus}$]':mass_i, 
      'Current Mass [M$_{\oplus}$]':mass_c,
      'Final Mass [M$_{\oplus}$]':mass_f,
      'Initial Radius [R$_{\oplus}$]':radius_i, 
      'Current Radius [R$_{\oplus}$]':radius_c,
      'Final Radius [R$_{\oplus}$]':radius_f,
      'Initial Core Fraction':core_i, 
      'Current Core Fraction':core_c,
      'Final Core Fraction':core_f,
      'Initial Envelope Fraction':env_i, 
      'Current Envelope Fraction':env_c,
      'Final Envelope Fraction':env_f,
      'Initial Envelope Composition':env_comp_i, 
      'Current Envelope Composition':env_comp_c,
      'Final Envelope Composition':env_comp_f,
      'Mass Lost [M$_{\oplus}$]':dmass
      })

df.to_pickle("./MC_save_secondrun.pkl")

