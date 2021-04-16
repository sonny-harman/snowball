import pickle as pkl
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import math
import re
import os

data = {}
runno = 0
runs = [f for f in os.listdir('saved_singleruns') if os.path.isfile(os.path.join('saved_singleruns',f))]

for runid in runs:
      with open('saved_singleruns/'+runid,'rb') as runfile:
            (m_planet, r_planet, a_planet, core_frac, core_den, albedo, m_star, age_star, l_star, 
                  d_star, age, envelope_species, c_f_t, e_f_t, r_p_t, m_p_t, mescape_flux_t, e_comp_t) = pkl.load(runfile)
      data[runid] = {'m_planet':m_planet, 'r_planet':r_planet, 'core_frac':core_frac, 
            'core_den':core_den, 'albedo':albedo, 'm_star':m_star, 'age_star':age_star, 
            'l_star':l_star, 'd_star':d_star, 'age':age, 'envelope_species':envelope_species, 
            'c_f_t':c_f_t, 'e_f_t':e_f_t, 'r_p_t':r_p_t, 'm_p_t':m_p_t, 
            'mescape_flux_t':mescape_flux_t, 'e_comp_t':e_comp_t}

fancy_species = [re.sub(r'(\d)','$_{\g<1>}$',species) for species in envelope_species]
colors = ['k','y','b','r']

fig,ax = plt.subplots(len(runs),1)
plt.subplots_adjust(wspace=0,hspace=0.1)
for run in data:
      single_run = data[run]
      x = single_run['age']
      y = single_run['e_comp_t']
      masses = {s:[y[i][envelope_species.index(s)] for i in range(len(y))] for s in envelope_species}
      for s in envelope_species:
            si = envelope_species.index(s)
            ax[runno].plot(x,masses[s],label=fancy_species[si],color=colors[si])
      ax[runno].plot([single_run['age_star'],single_run['age_star']],[0,1],ls='--',color='gray')
      ax[runno].set_yscale('symlog',linthresh=0.01)
      ax[runno].set_ylim(0,1)
      ax[runno].set_xlim(min(x),max(x))
      ax[runno].set_ylabel('Mass Mixing Ratio',fontsize=16)
      runno+=1

for i in range(len(runs)-1):
      ax[i].tick_params(axis='x',which='both',direction='in',labelbottom='off')
      ax[i].set_xticklabels([])
      for label in (ax[i].get_yticklabels()):
            label.set_fontsize(16)
for label in (ax[-1].get_yticklabels() + ax[-1].get_xticklabels()):
      label.set_fontsize(16)
ax[-1].set_xlabel('Time [Gyr]',fontsize=16)
ax[0].legend(loc='lower left',fontsize=16)
plt.show()
