import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from math import log10,floor
import corner
import copy
from planet import envelope_species

tickf = ticker.FuncFormatter(lambda y, _: '{:.2g}'.format(y))

fname = "MC_save_secondrun.pkl"
fnamefirst = fname.split('.')[0]
if not os.path.exists('./saved_plots/'+fnamefirst):
      os.mkdir('./saved_plots/'+fnamefirst)

df0 = pd.read_pickle("./"+fname)

#still hunting unphysical mass-radius combinations
df = df0[df0["Initial Radius [R$_{\oplus}$]"] > 1.55]
df = df[df["Current Radius [R$_{\oplus}$]"] > 1.55]
df = df[df["Final Radius [R$_{\oplus}$]"] > 1.55]

df.info(verbose=True)

#df.plot(kind='scatter',x='Current Mass [M$_{\oplus}$]',xlabel='Present-day Mass [M$_{\oplus}$]',y='Current Envelope Fraction',ylabel='Present-Day Envelope Fraction')
#plt.show()

dim1,dim2 = df.shape
#Envelope composition is in terms of envelope_species 
# 2   Loss Efficiency                1000 non-null   object 
# 15  Initial Envelope Composition   1000 non-null   object 
# 16  Current Envelope Composition   1000 non-null   object 
# 17  Final Envelope Composition     1000 non-null   object 
leff = df['Loss Efficiency']
Heff = [eff[0] for eff in leff]
iec = df['Initial Envelope Composition']
cec = df['Current Envelope Composition']
fec = df['Final Envelope Composition']
for s in envelope_species:
      df[f"Initial {s} Mass Fraction"] = [x[envelope_species.index(s)] for x in iec]
      df[f"Current {s} Mass Fraction"] = [x[envelope_species.index(s)] for x in cec]
      df[f"Final {s} Mass Fraction"] = [x[envelope_species.index(s)] for x in fec]

df['H2 Loss Efficiency'] = Heff
#for s in ['H2', 'H2O', 'CO2']:
#      df[f'Delta {s} [M$_{\oplus}$]'] = (df['Current {s} Mass Fraction'].multiply() ).subtract(df['Initial {s} Mass Fraction'].multiply())

df = df.drop(columns=['Loss Efficiency','Initial Envelope Composition','Current Envelope Composition','Final Envelope Composition','Initial O Mass Fraction'])
df.info(verbose=True)

df2 = copy.deepcopy(df)
dropem = []
limited_list = ['Initial Core Fraction', 'Initial Envelope Fraction', 'Initial Mass [M$_{\\oplus}$]', 'Initial H2O Mass Fraction', 'Initial CO2 Mass Fraction', 'H2 Escape Efficiency', 'Current H2 Mass Fraction', 'Current H2O Mass Fraction', 'Current O Mass Fraction', 'Mass Lost [M$_{\\oplus}$]']
#for key in df2:
#      if 'Final' in key or 'He Mass' in key or 'Core' in key:
#            dropem.append(key)
#for d in dropem:
#      df2.pop(d)
labels = df2.keys()
for l in labels:
      if l not in limited_list:
            df2.pop(l)
labels = df2.keys()
print(df2.mean(axis=0))
print(df2.var(axis=0))
dscrp = df2.describe()
print(dscrp.transpose())
#print('Covariance matrix:')
#print(df2.cov())
#print('Correlation matrix:')
#print(df2.corr())
corrvals = df2.corr()
small = []; medium = []; large = []
for x in df2:
      for y in df2:
            val = corrvals[x][y]
            if x != y and val != 'NaN':
                  if 0.1 <= abs(val) < 0.3:
                        small.append([x,y,val])
                        corrvals[x][y] = 'NaN'
                        corrvals[y][x] = 'NaN'
                  elif 0.3 <= abs(val) < 0.5:
                        medium.append([x,y,val])
                        corrvals[x][y] = 'NaN'
                        corrvals[y][x] = 'NaN'
                  elif 0.5 <= abs(val):
                        large.append([x,y,val])
                        corrvals[x][y] = 'NaN'
                        corrvals[y][x] = 'NaN'
      #print(corrvals[x])
#print('small = ',small)
#print('medium = ',medium)
#print('large = ',large)
#for pair in large:
      #print(pair)
      #if 'Initial' in pair[0] or 'Initial' in pair[1]:
      #      print(pair)

numpyarr = df2.to_numpy()
#figure = corner.corner(numpyarr, labels=labels)
#plt.show()
#figure.savefig('saved_plots/corner_MC_v2.pdf',dpi=300)
#exit()

#exclude_vars = 

yvar = 'Current H2 Mass Fraction'
#print('Testing where the non-zero oxygen cases live; yvar = '+yvar)
#df = df[df[yvar] <= 0.] #only look at cases with non-zero oxygen amounts

#df.plot(kind='scatter',x='Current Mass [M$_{\oplus}$]',xlabel=r'Present-day Mass [M$_{\oplus}$]',y='Current H2O Mass Fraction',ylabel=r'Present-Day H$_{2}$O Mass Fraction')
#plt.yscale('log')
#plt.show()
mmax = max(df['Initial Envelope Fraction'])
mmin = min(df['Initial Envelope Fraction'])
def_msize = 0.01
#sizing = lambda s: def_msize*(1+((log10(100*s)-log10(100*mmin))/(log10(100*mmax)-log10(100*mmin)))*10000)
sizing = lambda s,min,max: def_msize*(1+((100*s)-(100*min))/((100*mmax)-(100*min))*5000)
size = [sizing(s,mmin,mmax) for s in df['Current Mass [M$_{\oplus}$]']]#df['Initial Envelope Fraction']]
print('Max H2O = ',max(df['Current H2O Mass Fraction']))

amtH2 = [i*ie*im-c*e*m for c,e,m,i,ie,im in zip(df['Current H2 Mass Fraction'],df['Current Envelope Fraction'],df['Current Mass [M$_{\oplus}$]'],df['Initial H2 Mass Fraction'],df['Initial Envelope Fraction'],df['Initial Mass [M$_{\oplus}$]'])]
initH2 = [i*ie*im for c,e,m,i,ie,im in zip(df['Current H2 Mass Fraction'],df['Current Envelope Fraction'],df['Current Mass [M$_{\oplus}$]'],df['Initial H2 Mass Fraction'],df['Initial Envelope Fraction'],df['Initial Mass [M$_{\oplus}$]'])]
presH2 = [c*e*m for c,e,m,i,ie,im in zip(df['Current H2 Mass Fraction'],df['Current Envelope Fraction'],df['Current Mass [M$_{\oplus}$]'],df['Initial H2 Mass Fraction'],df['Initial Envelope Fraction'],df['Initial Mass [M$_{\oplus}$]'])]

anyH2O = ['k' if s*e*m>0.002 else 'r' for s,e,m in zip(df['Current H2O Mass Fraction'],df['Current Envelope Fraction'],df['Current Mass [M$_{\oplus}$]'])]
lt_one_ocean = ['k' if s*e*m>0.002 else 'r' for s,e,m in zip(df['Current H2O Mass Fraction'],df['Current Envelope Fraction'],df['Current Mass [M$_{\oplus}$]'])]
CO2plusH2 = ['k' if (c>0.0005 and c<0.3 and h>0.00999) else 'r' for c,h in zip(df['Current CO2 Mass Fraction'],df['Current H2 Mass Fraction'])]
mu_bar = [1/((2*h + 4*he + 18*h2o + 44*co2 + 16*o)-1.E-10) for h,he,h2o,co2,o in zip(df['Current H2 Mass Fraction'],df['Current He Mass Fraction'],df['Current H2O Mass Fraction'],df['Current CO2 Mass Fraction'],df['Current O Mass Fraction'])]
CO2 = [c for c in df['Current CO2 Mass Fraction']]
co2_vmr = [CO2[i]*mu_bar[i]/44. for i in range(len(mu_bar))]
someCO2 = ['k' if (s < 0.3 and s > 0.0005) else 'r' for s in co2_vmr]
anyCO2 = [c for c in co2_vmr if c > 1E-6] #df[df['Current CO2 Mass Fraction'] > 1E-6]['Current CO2 Mass Fraction']
sCO2 = [c for c in co2_vmr if c < 0.3] #df[df['Current CO2 Mass Fraction'] < 0.3]['Current CO2 Mass Fraction']
sCO2 = [c for c in sCO2 if c > 0.0005] #sCO2[sCO2 > 0.0005]
print(f"    There are {len(sCO2)} cases where there is detectable CO2, and {len(anyCO2)} cases where there is > 1ppm CO2.")
anyO = df[df['Current O Mass Fraction'] == 0]['Current O Mass Fraction']
print(f"    There are {CO2plusH2.count('k')} cases where there is CO2 and enough H2 to see it, in 10k cases.")
print(f"    There are {anyH2O.count('r')} cases where there is less than 1 ocean of water, and {len(anyO)} cases where there is no oxygen.")
#size = [s*def_msize*300 for s in df['Initial Envelope Fraction']]
#ax = df.plot(kind='scatter',x='Current Mass [M$_{\oplus}$]',xlabel=r'Present-day Mass [M$_{\oplus}$]',y='Current O Mass Fraction',ylabel=r'Present-Day O Mass Fraction',c='Initial H2O Mass Fraction',cmap='Blues',edgecolor=anyH2O,linewidth=1,s=size,logy=True,fontsize=16)
print(f"Earth masses of H2 lost: {min(amtH2):.2e}-{np.mean([x for x in amtH2 if x>0]):.2e}-{max(amtH2):.2e}")

#df.to_excel('output.xlsx')

#negH2 = [amtH2.index(x) for x in amtH2 if x<0]
spec_cases = [amtH2.index(x) for x in amtH2 if x<0]
print(df.mean(axis=0))
print(df.var(axis=0))
dscrp = df.describe()
print(dscrp.transpose())
#for x in spec_cases:
#      print(df.loc[x,:])
#      print(df['Initial H2 Mass Fraction'][x],'   ',df['Initial Envelope Fraction'][x],df['Initial Mass [M$_{\oplus}$]'][x],df['Current H2 Mass Fraction'][x],df['Current Envelope Fraction'][x],df['Current Mass [M$_{\oplus}$]'][x])
#print(df['Initial H2 Mass Fraction'][negH2])
badH2 = ['r' if c>i else 'b' for i,c in zip(df['Initial H2 Mass Fraction'],df['Current H2 Mass Fraction'])]
print(f"There are {badH2.count('r')} cases where the current H2 mass fraction is larger than the initial.")

#fig = plt.figure()
##plt.scatter(initH2,presH2,c=badH2)
#perc_env_change = (df['Initial Envelope Fraction'] - df['Current Envelope Fraction'])/df['Initial Envelope Fraction'] 
#print("Mean of envelope mass fraction in % change = ",np.mean(perc_env_change))
#plt.scatter(df['Current Mass [M$_{\oplus}$]'],perc_env_change,c=df['Initial Envelope Fraction'])
#plt.xscale('log')
#plt.yscale('log')
#plt.scatter(df['Initial H2 Mass Fraction'],df['Current H2 Mass Fraction'],c=badH2)
#plt.plot(df['Initial Envelope Fraction'],df['Current Envelope Fraction'],'o')
#plt.plot([0,1],[0,1],'k')
#plt.plot(df['Current H2O Mass Fraction'],amtH2,'o')
#plt.xlabel('Initial H2 Mass Fraction')
#plt.ylabel('Current H2 Mass Fraction')
#plt.xlabel('Initial H2 Planet Mass Fraction')
#plt.ylabel('Current H2 Planet Mass Fraction')
#plt.xlabel('Current Mass')
#plt.ylabel('% Change in Envelope Fraction')
#plt.xlabel('Initial Envelope Fraction')
#plt.ylabel('Current Envelope Fraction')
#plt.show()

#exit()

#df.plot(kind='scatter',x='Initial Radius [R$_{\oplus}$]',xlabel='Initial Radius',y='Current Radius [R$_{\oplus}$]',ylabel='Current Radius')
#df.plot(kind='scatter',x='Initial Mass [M$_{\oplus}$]',xlabel='Initial Mass',y='Initial Radius [R$_{\oplus}$]',ylabel='Initial Radius')
#plt.show()

#yvar = 'Current Mass [M$_{\oplus}$]'
#limy = [df[yvar].min(), df[yvar].max()] 
limy = [3.E-5, 1]#[df[df[yvar] > 1.E-9][yvar].min(),1]
print('y limits = ',limy)
quants = [0.01, 0.1, 0.33, 0.5, 0.68, 0.9, 0.99]
yvar_q = [df.quantile(x) for x in quants]
yvar_quantile = [y[yvar] for y in yvar_q]
yvar_quant_lookup = {}
for p in quants[::-1]:
      if yvar_quantile[quants.index(p)] not in yvar_quant_lookup.values():
            yvar_quant_lookup[p] = yvar_quantile[quants.index(p)]
print(yvar_quant_lookup)

H2noH2O = len(df[(df['Current H2 Mass Fraction'] > 0) & (df['Current H2O Mass Fraction'] <= 0)])/len(df)
H2andH2O = len(df[(df['Current H2 Mass Fraction'] > 0) & (df['Current H2O Mass Fraction'] > 0)])/len(df)
H2overH2O = len(df[(df['Current H2 Mass Fraction'] > 0) & (df['Current H2O Mass Fraction'] > 0)])/len(df[df['Current H2O Mass Fraction'] > 0])
H2OoverH2 = len(df[(df['Current H2 Mass Fraction'] > 0) & (df['Current H2O Mass Fraction'] > 0)])/len(df[df['Current H2 Mass Fraction'] > 0])
print('Fraction of runs where there is H2 and H2O: ',H2andH2O)
print('Fraction of runs with H2O where there is H2: ',H2overH2O)
print('Fraction of runs with H2 where there is H2O: ',H2OoverH2)
print('Fraction of runs where there is H2 and no H2O: ',H2noH2O)

def plot_it(x,y,s,c,xl,yl,sl,cl,edge='k',map='Blues',log=True,scale=1,shift=0,units=1):
      fig,ax = plt.subplots(figsize=(6.5,4.5),dpi=100)
      size = [sizing(si,min(s),max(s))*scale+shift for si in s]
      p1 = plt.scatter(x,y,c=c,s=size,edgecolor=edge,cmap=map)
      #add a size key
      m1 = plt.scatter([],[],s=min(size),color='gray')
      m2 = plt.scatter([],[],s=np.mean(size),color='gray')
      m3 = plt.scatter([],[],s=max(size),color='gray')
      legend_markers = [m1, m2, m3]
      labels = [str(round(x*units,1)) for x in [min(s),np.mean(s),max(s)]]
      ax.legend(handles=legend_markers,labels=labels,loc='best',scatterpoints=1,fontsize=16)
      #add colorbar
      #fig.colorbar(ax,ax=ax,cmap='Blues')
      cbar_axes = plt.colorbar(p1,ax=ax)
      cbar_axes.ax.set_ylabel(cl,rotation=90)
      cbar_axes = ax.figure.axes[-1].yaxis.label.set_size(16)
      cbar_axes = ax.figure.axes[-1].tick_params(labelsize=16)
      if log==True:
            if min(x)/max(x) < 0.05:
                  plt.xscale('log')
            if min(y)/max(y) < 0.05:
                  plt.yscale('log')
      plt.xlim(min(x),max(x))
      #ax.axes.
      plt.ylim(min(y)+1E-5,max(y))
      plt.xlabel(xl,fontsize=16)
      plt.ylabel(yl,fontsize=16)
      plt.tight_layout()
      plt.savefig('./saved_plots/plot_it/'+xl+'_v_'+yl+'_by_'+cl+'.pdf')
      plt.show()

zvar = 'Initial Envelope Fraction'
size = [sizing(s,min(df[zvar]),max(df[zvar])) for s in df[zvar]]
#plot_it(df['Current Mass [M$_{\oplus}$]'],df['Current O Mass Fraction'],df[zvar],df['Initial H2O Mass Fraction'],
#      r'Current Mass [M$_{\oplus}$]','Current O Envelope Mass Fraction',zvar,r'Initial Envelope H$_{2}$O Mass Fraction',
#      edge=lt_one_ocean,map='Blues',scale=2,shift=20,units=100)
plot_it(Heff,df['Current O Mass Fraction'],df[zvar],df['Initial H2O Mass Fraction'],
      r'H$_{2}$ Escape Efficiency','Current O Envelope Mass Fraction',zvar,r'Initial Envelope H$_{2}$O Mass Fraction',
      edge=lt_one_ocean,map='Blues',log=False,scale=2,shift=20,units=100)

#plot the probability density functions of present-day oxygen, water, and hydrogen amounts as a function of the following variables:
# H2 escape efficiency, PD CO2/O2/H2O/H2, 
print('Covariance matrix:')
#print(df.cov())
print('Correlation matrix:')
#print(df.corr())
corrvals = df.corr()
small = []; medium = []; large = []
for x in df:
      for y in df:
            val = corrvals[x][y]
            if x != y and val != 'NaN':
                  if 0.1 <= abs(val) < 0.3:
                        small.append([x,y,val])
                        corrvals[x][y] = 'NaN'
                        corrvals[y][x] = 'NaN'
                  elif 0.3 <= abs(val) < 0.5:
                        medium.append([x,y,val])
                        corrvals[x][y] = 'NaN'
                        corrvals[y][x] = 'NaN'
                  elif 0.5 <= abs(val):
                        large.append([x,y,val])
                        corrvals[x][y] = 'NaN'
                        corrvals[y][x] = 'NaN'
      #print(corrvals[x])
#print('small = ',small)
#print('medium = ',medium)
#print('large = ',large)
for pair in large:
      if 'Initial' in pair[0] or 'Initial' in pair[1]:
            print(pair)


#X/Y/size/H2O plot
#fig1 = plt.figure()
xvar = 'Current Mass [M$_{\oplus}$]'
yvar = 'Current O Mass Fraction'
ax = df.plot(kind='scatter',x=xvar,xlabel=xvar,y=yvar,ylabel=yvar,c='Initial H2O Mass Fraction',cmap='Blues',edgecolor=anyH2O,linewidth=1,s=size,fontsize=16,figsize=(10,7))
#add a size key
#Invalid value encountered in the sqrt - negative number in min term?
mindat = min(df['Current Mass [M$_{\oplus}$]'])
maxdat = max(df['Current Mass [M$_{\oplus}$]'])
m1 = plt.scatter([],[],s=sizing(0.5,mindat,maxdat),color='gray')#,s=sizing(0.5,max(0,min(df['Current Mass [M$_{\oplus}$]'])),max(df['Current Mass [M$_{\oplus}$]'])),color='gray')
m2 = plt.scatter([],[])#,s=sizing(1.9,min(df['Current Mass [M$_{\oplus}$]']),max(df['Current Mass [M$_{\oplus}$]'])),color='gray')
m3 = plt.scatter([],[])#,s=sizing(6.4,min(df['Current Mass [M$_{\oplus}$]']),max(df['Current Mass [M$_{\oplus}$]'])),color='gray')
legend_markers = [m1, m2, m3]
labels = [str(round(x,-int(floor(log10(abs(x))))))+'%' for x in [0.01,1,100]]
ax.legend(handles=legend_markers,labels=labels,loc='best',scatterpoints=1,fontsize=16)
#add colorbar
cbar_axes = ax.figure.axes[-1].yaxis.label.set_size(16)
cbar_axes = ax.figure.axes[-1].tick_params(labelsize=16)
plt.yscale('log')
plt.ylim([1e-4,1])
plt.tight_layout()
plt.show()
exit()

#varlist = [v if v not in exclude_vars for v in df]
for var in df:
      limx = [max(1.E-9,df[df[var] > 0.][var].min()), df[var].max()]
      #ax = df.plot(kind='scatter',x=var,xlabel=var,y=yvar,ylabel=yvar,c='Initial H2O Mass Fraction',cmap='Blues',edgecolor=anyH2O,linewidth=1,s=size,fontsize=16,figsize=(10,7))
      #ax = df.plot(kind='scatter',x=var,xlabel=var,y=yvar,ylabel=yvar,c='Current Mass [M$_{\oplus}$]',cmap='Blues',edgecolor=CO2plusH2,linewidth=1,s=size,fontsize=16,figsize=(10,7))
      ax = df.plot(kind='scatter',x=var,xlabel=var,y=yvar,ylabel=yvar,c='Current CO2 Mass Fraction',cmap='Blues',edgecolor=CO2plusH2,linewidth=1.5,s=size,fontsize=16,figsize=(10,7))
      plt.yscale('symlog',linthresh=limy[0])
      if df[var].min() <= 0:
            plt.xscale('symlog',linthresh=limx[0])
            print(var,f'xrange = {[df[var].min(),df[var].max()]}, linear threshold = {limx[0]}')
      else:
            plt.xscale('log')
            print(var,f'xrange = {[df[var].min(),df[var].max()]}, logscale')
      
      #add a size key
      m1 = plt.scatter([],[],s=sizing(0.5,min(df['Current Mass [M$_{\oplus}$]']),max(df['Current Mass [M$_{\oplus}$]'])),color='gray')
      m2 = plt.scatter([],[],s=sizing(1.9,min(df['Current Mass [M$_{\oplus}$]']),max(df['Current Mass [M$_{\oplus}$]'])),color='gray')
      m3 = plt.scatter([],[],s=sizing(6.4,min(df['Current Mass [M$_{\oplus}$]']),max(df['Current Mass [M$_{\oplus}$]'])),color='gray')
      #m1 = plt.scatter([],[],s=sizing(0.01),color='gray')
      #m2 = plt.scatter([],[],s=sizing(1),color='gray')
      #m3 = plt.scatter([],[],s=sizing(100),color='gray')
      #labels = [str(round(x,-int(floor(log10(abs(x))))))+'%' for x in [0.01,1,100]]
      legend_markers = [m1, m2, m3]
      labels = [str(x) for x in [0.5,1.9,6.4]]
      ax.legend(handles=legend_markers,labels=labels,loc='best',scatterpoints=1,fontsize=16)
      
      #add colorbar
      cbar_axes = ax.figure.axes[-1].yaxis.label.set_size(16)
      cbar_axes = ax.figure.axes[-1].tick_params(labelsize=16)
      
      #point out where the axis transitions from log to linear
      ax.text(max(0,df[var].min()), limy[0], f"$\\approx$ ", fontsize=16, horizontalalignment='center', verticalalignment='center')
      ax.text(max(1.E-9,df[df[var] > 0.][var].min()), 0, f"$\\approx$ ", rotation=90,fontsize=16, horizontalalignment='center', verticalalignment='center')
      
      #change plotting preferences
      #ax.xaxis.label.set_size(16)
      #ax.yaxis.label.set_size(16)
      plt.ylim([0,1])
      plt.xticks(fontsize=16,fontweight='bold')
      plt.yticks(fontsize=16,fontweight='bold')
      plt.xlim([max(0.,df[var].min()),df[var].max()])
      plt.xlabel(var,fontsize=16,fontweight='bold')
      plt.ylabel(yvar,fontsize=16,fontweight='bold')

      
      plt.tight_layout()
      plt.show()
      #plt.savefig('./saved_plots/'+fnamefirst+'/'+yvar.replace(' ','_')+'_vs_'+var.split('[')[0]+ax.get_ylabel().replace(' ','_')+fnamefirst+'.png',bbox_inches='tight',dpi=300)
      plt.close()



exit() 
#Grab just the ~0.1 H escape efficiency values
#downselect = df[(df['H2 Loss Efficiency'] > 0.35) & (df['H2 Loss Efficiency'] < 0.5)]
#downselect = downselect.sort_values(by=['Current Mass [M$_{\oplus}$]'])
##print(min(downselect['H2 Loss Efficiency']),np.mean(downselect['H2 Loss Efficiency']),max(downselect['H2 Loss Efficiency']))
#xvar = downselect['Current Mass [M$_{\oplus}$]']
#y = downselect['Current O Mass Fraction']#.rolling(window=10).mean()
##yvar = [sum(y[max(0,i-2):min(len(y),i+2)])/5 for i in range(len(y))]
#yvar = y.rolling(20).mean()#y.ewm(span=10,adjust=False).mean()
#ax.plot(xvar,yvar)
#plt.yscale('log')
##cbar = plt.colorbar()
##cbar.set_label('Initial H$_{2}$O Mass Fraction')

#pd.plotting.scatter_matrix(df,alpha=0.5)
#plt.show()

df.plot(kind='scatter',x='Initial H2O Mass Fraction',xlabel=r'Initial H$_{2}$O Mass Fraction',y='Current O Mass Fraction',ylabel=r'Present-Day O Mass Fraction',s=40,c=Heff,cmap='viridis')
plt.plot([0,1],[0,16/18])
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-4,1])
plt.xlim([1.E-4,1])
plt.tight_layout()
plt.show()


#x = df['Current Mass [M$_{\oplus}$]']
#y = df['Current O Mass Fraction'] 
#z = df['Initial H2O Mass Fraction']
#z = Heff

#dfnew = pd.DataFrame({'Mass':x,'O':y,'H2 loss efficiency':z})
#df.plot(kind='hexbin',x='Current Mass [M$_{\oplus}$]',xlabel=r'Present-day Mass [M$_{\oplus}$]',y='Current O Mass Fraction',ylabel=r'Present-Day O Mass Fraction',C='Initial H2O Mass Fraction',reduce_C_function=np.mean,gridsize=40)
#dfnew.plot(kind='hexbin',x='Mass',xlabel=r'Present-day Mass [M$_{\oplus}$]',y='O',ylabel=r'Present-Day O Mass Fraction',C='H2 loss efficiency',reduce_C_function=np.mean,gridsize=25)
#plt.show()


#df.cov #yields the pairwise covariance of columns
#df.plot.density()
#plt.show()
