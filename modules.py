#This file contains modules for the hydrodynamic escape code to call
#     when necessary; users are encouraged to think before modifying
#     the underlying physics of the model.

#mubar is the befault mean molecular weight of the atmosphere
from planet import mubar
from constants import kb,Rconst,m_H
from math import log10
from bisect import bisect_left,bisect_right
from os import linesep,listdir

#Handy lambdas for finding the nearest,next above, and next below value in a list.
find_nearest = lambda vector,value : min(range(len(vector)), key = lambda i: abs(vector[i]-value))
find_above = lambda vector,value : bisect_right(vector,value)
find_below = lambda vector,value : max(bisect_left(vector,value)-1,0)

#lambda function for finding the luminosity (in solar) for a 0.2<M<0.85 M_solar star (Cuntz & Wang, 2018)
f_lum_CW18 = lambda M : M**(-141.7*M**4. + 232.4*M**3. - 129.1*M**2. + 33.29*M + 0.215)

def generic_diffusion(minor,major,T):
      #generic diffusion of a minor species with mean molecular weight 
      #     *minor* (amu) against a major species with mmw *major* at 
      #     temperature T (K); from Banks and Kockarts (1973)/ Kasting
      #     units converted from /cm/s to /m/s
      b = 1.52E20*(1./minor + 1./major)**0.5*(T**0.5)
      bLB15 = 4.8E19*T**0.75 #/m/s, from Luger & Barnes (2015) and Zahnle & Kasting (1986)
      return b

def crossover_mass(T,flux,gravity,X,mass=1.,mu=mubar):
      b=generic_diffusion(mass,mu,T)
      #crossover_mass(T_eq_t[-1],eflux,grav_t[-1],2.*e_comp_t[-1][0],mass=1.,mu=mu_t[-1])
      #calculate the crossover mass at a given escape flux based on
      #     the temperature T (K), flux (molecules/m2/s), gravity 
      #     (m/s2), mixing ratio of escaping species X, and binary 
      #     diffusion parameterizations given by generic_diffusion 
      #     (/m/s).
      cm = mass + kb*T*flux/(b*m_H*gravity*X)
      return cm

def read_baraffe(mass):
      #This loads the temperature and luminosity evolution
      #     for a star with the given *mass* solar masses
      #     from Baraffe et al. (2015), online at:
      #     http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/
      masses = []; ages = {}; teff = {}; lums = {}; gravs = {}
      stop_after_this = False
      with open('data/BHAC15_tracks+structure.txt','r') as fin:
            lines_after_header = fin.readlines()[45:]
      #Variables in Baraffe et al. file:
      #! M/Ms    log t(yr)  Teff     L/Ls    g   R/Rs Log(Li/Li0) log Tc  log ROc   Mrad     Rrad       k2conv      k2rad
      #print(lines_after_header[0])
      for line in lines_after_header:
            if line[0] == '!' or line[0] in linesep:
                  #only read the data until there are values on both sides of the input mass
                  if stop_after_this == True:
                        break
                  #reset current mass counter between blocks
                  current_mass = 0.
                  #ignore header lines for each block
                  continue
            else:
                  l = line.split()
                  m = float(l[0]); a = 10**float(l[1])/1.E9; t = float(l[2]); lum = 10**float(l[3]); g = float(l[4])
                  key = round(m,4)
                  if current_mass == 0.:
                        current_mass = key
                        if current_mass >= mass:
                              m1 = masses[-1]; m2 = current_mass
                              stop_after_this = True
                        ages[key] = [a]; teff[key] = [t]
                        lums[key] = [lum]; gravs[key] = [g]
                        masses.append(key)
                  elif current_mass == m:
                        ages[key].append(a); teff[key].append(t)
                        lums[key].append(lum); gravs[key].append(g)
      age1 = ages[m1]; age2 = ages[m2]
      lum1 = lums[m1]; lum2 = lums[m2]
      return age1,age2,lum1,lum2,m1,m2

def read_baraffe_grid():
      #This loads the temperature and luminosity evolution
      #     for all stars from Baraffe et al. (2015), online at:
      #     http://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/
      mass = []; age = []; temp = []; luminosity = []; grav = []
      count = 0 # just to make sure there's a count of the number of stellar masses
      with open('data/BHAC15_tracks+structure.txt','r') as fin:
            lines_after_header = fin.readlines()[45:]
      for line in lines_after_header:
            if line[0] == '!' or line[0] in linesep:
                  continue
            else:
                  l = line.split() #currently ignoring other data
                  m = float(l[0]); a = 10**float(l[1])/1.E9; t = float(l[2]); lum = 10**float(l[3]); g = float(l[4])
                  if m not in mass:
                        count += 1
                  mass.append(m); age.append(a); temp.append(t); luminosity.append(lum); grav.append(g)
      print(f"Number of stars in grid: {count:2.0f}")
      return mass,age,temp,luminosity,grav

def read_hidalgo(stellar_mass):
      #This loads the temperature and luminosity evolution 
      #of bookending stars from Hidalgo et al. (2018), online at:
      # http://basti-iac.oa-abruzzo.inaf.it
      #http://basti-iac.oa-abruzzo.inaf.it/cgi-bin/track-get.py?alpha=P00&grid=P00O1D1E1Y247&metal=FEHp006&imetal=&imetalh=&mass=all&imass=&bcsel=HR
      #[α/Fe]=+0.0, Z=0.017210, He Abundance=0.26948852, Oversh.=Yes, Diffusion=Yes, Mass loss=0.3
      datafiles = [f for f in listdir('data/BaSTI_stars/')]
      masses = [float(n[:2]+'.'+n[2:4]) for n in datafiles]
      masses.sort()
      f1 = [masses[find_below(masses,stellar_mass)], masses[find_above(masses,stellar_mass)]]
      f2 = [f"{f:0>5.2f}".replace('.','') for f in f1]
      for file_start in f2:
            match = [i for i in datafiles if file_start in i]
            with open('data/BaSTI_stars/'+match[0],'r') as fin:
                  age = []; luminosity = []
                  lines = fin.readlines()[5:]
                  for line in lines:      
                        l = line.split() #currently ignoring other data
                        a = 10**float(l[0])/1.E9; lum = 10**float(l[1]); t = float(l[2])
                        age.append(a); luminosity.append(lum)
            if f2.index(file_start) == 0:
                  a1 = age; l1 = luminosity
            else:
                  a2 = age; l2 = luminosity
      return a1,a2,l1,l2,f1[0],f1[1]

def read_thermo(name):
      coeffs = [[],[]]
      with open('data/thermo.dat','r') as fin:
            lines_after_header = fin.readlines()[10:]
      for line in lines_after_header:
            lineno = lines_after_header.index(line)
            l = line.split()
            if not l:
                  continue
            elif l[0] == name:
                  print('Found '+name+': '+l[0])
                  for i in range(2):
                        line1 = lines_after_header[lineno+1+i*2].split()
                        line2 = lines_after_header[lineno+2+i*2].split()
                        coeffs[i] = [float(var)*Rconst for var in (line1 + line2)]
      if lineno == len(lines_after_header):
            print(name+' not found in read_thermo!')
      return coeffs

def read_pt_profile():
      p = []; alt = []; t = []
      with open('data/clima_last_0.01H2_0.99H2O_200KStrat.dat','r') as fin:
            lines_after_header = fin.readlines()[1:]
      for l in lines_after_header:
            p.append(float(l.split()[1])*1E5)
            alt.append(float(l.split()[2])*1e3)
            t.append(float(l.split()[3]))
      return p,alt,t

def new_func():
      return
