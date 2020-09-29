from math import log10
from numpy import logspace,linspace
from main import run_escape
import pandas as pd
import numpy as np
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
ntimes = 10

MC_mass = logspace(log10(0.6),log10(6.4),ntimes,endpoint=True)
MC_age = [7.9]*ntimes
MC_lum = [0.02629]*ntimes
MC_radius = [1.673]*ntimes
MC_envelope = [[0.1,0.01,0.88,0.0,0.01]]*ntimes
MC_eff = [[0.1,0.01,0.01,0.01,0.01]]*ntimes

age = []; lum=[]; mass = []; radius = []; envelope = []; eff = []

for i in range(ntimes):
      MC_params = {'age':MC_age[i],
                  'lum':MC_lum[i],
                  'mass':MC_mass[i],
                  'radius':MC_radius[i],
                  'envelope':MC_envelope[i],
                  'efficiency':MC_eff[i]
                  }
      
      mi,mc,mf,ri,rc,rf,ci,cc,cf,ei,ec,ef,eci,ecc,ecf,mass_change = run_escape(
                        MC_age[i],MC_lum[i],MC_mass[i],MC_radius[i],MC_envelope[i],MC_eff[i])
      
      #df = pd.DataFrame({'Current Age [Gyr]':MC_age[i],
      #                  'Luminosity [L/Lsol]':MC_lum[i],
      #                  ''
      #                  'Initial Mass':mi, 
      #                  'Current Mass':mc,
      #                  ''
