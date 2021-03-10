import pickle as pkl
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

runs = ['2021-03-10 09_43_51.070318']
data = {}

for runid in runs:
      with open('saved_singleruns/'+runid+'.run','rb') as runfile:
            (m_planet, r_planet, a_planet, core_frac, core_den, albedo, m_star, age_star, l_star, 
                  d_star, age, envelope_species, c_f_t, e_f_t, r_p_t, m_p_t, mescape_flux_t, e_comp_t) = pkl.load(runfile)
      
#do stuff
