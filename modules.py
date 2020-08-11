#This file contains modules for the hydrodynamic escape code to call
#     when necessary; users are encouraged to think before modifying
#     the underlying physics of the model.

#mubar is the befault mean molecular weight of the atmosphere
from planet import mubar
from constants import kb,Rconst,m_H
from math import log10
from bisect import bisect_left,bisect_right
from os import linesep

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
      with open('data/BHAC15_tracks+structure.txt','r') as fin:
            lines_after_header = fin.readlines()[45:]
      #Variables in Baraffe et al. file:
      #! M/Ms    log t(yr)  Teff     L/Ls    g   R/Rs Log(Li/Li0) log Tc  log ROc   Mrad     Rrad       k2conv      k2rad
      #print(lines_after_header[0])
      for line in lines_after_header:
            if line[0] == '!' or line[0] in linesep:
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
                        ages[key] = [a]; teff[key] = [t]
                        lums[key] = [lum]; gravs[key] = [g]
                        masses.append(key)
                  elif current_mass == m:
                        ages[key].append(a); teff[key].append(t)
                        lums[key].append(lum); gravs[key].append(g)
      #print(masses)
      match = masses[min(range(len(masses)), key = lambda i: abs(masses[i]-mass+1.E-10))]
      message = "read_baraffe is returning values for a "+str(match)+" solar mass star."
      message+= "\n     - Age is in Gyr."
      message+= "\n     - Luminosity is L/L_sun"
      subage = ages[match]; sublum = lums[match]
      return match,subage,sublum,message

def read_baraffe_grid():
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
      limits = [min(age),max(age),min(mass),max(mass)]
      print(f"Number of stars in grid: {count:2.0f}")
      return mass,age,temp,luminosity,grav,limits

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

#Handy lambdas for finding the nearest,next above, and next below value in a list.
find_nearest = lambda vector,value : min(range(len(vector)), key = lambda i: abs(vector[i]-value))
find_above = lambda vector,value : bisect_right(vector,value)
find_below = lambda vector,value : max(bisect_left(vector,value)-1,0)

#lambda function for finding the luminosity (in solar) for a 0.2<M<0.85 M_solar star (Cuntz & Wang, 2018)
f_lum_CW18 = lambda M : M**(-141.7*M**4. + 232.4*M**3. - 129.1*M**2. + 33.29*M + 0.215)

def new_func():
      return
