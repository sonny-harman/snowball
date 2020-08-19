#This file contains modules for the hydrodynamic escape code to call
#     when necessary; users are encouraged to think before modifying
#     the underlying physics of the model.

#mubar is the befault mean molecular weight of the atmosphere
from planet import a_planet,hnu,envelope_species
from constants import kb,Rconst,m_H,h,c,G,r_Earth
from math import log10,exp
from bisect import bisect_left,bisect_right
from os import linesep,listdir
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.special import lambertw

#Handy lambdas for finding the nearest,next above, and next below value in a list.
find_nearest = lambda vector,value : min(range(len(vector)), key = lambda i: abs(vector[i]-value))
find_above = lambda vector,value : bisect_right(vector,value)
find_below = lambda vector,value : max(bisect_left(vector,value)-1,0)

#Noack et al. (2016) water/iron/silicate mass-radius relationships
AN = lambda xfe,xh2o: 7121 - 20.21*xfe + xh2o*(15.23 + 0.239*xfe)
ANFE = lambda xfe: 7121 - 20.21*xfe + (1-xfe)*(15.23 + 0.239*xfe)
ANH2O = lambda xh2o: 7121 + xh2o*(15.23)
CN = lambda xh2o: 0.2645 + 0.00048*xh2o
f_RpN = lambda mp,xfe,xh2o: 1000*AN(xfe,xh2o)*mp**CN(xh2o)/r_Earth

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

def crossover_mass(T,flux,gravity,X,mu,mass=1.):
      b=generic_diffusion(mass,mu,T)
      #calculate the crossover mass at a given escape flux based on
      #     the temperature T (K), flux (molecules/m2/s), gravity 
      #     (m/s2), mixing ratio of escaping species X, and binary 
      #     diffusion parameterizations given by generic_diffusion 
      #     (/m/s).
      cm = mass + kb*T*flux/(b*m_H*gravity*X)
      return cm

def planet_radius(mass,xfe,xh2o):
      radius = f_RpN(mass,xfe,xh2o)
      return radius

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

def read_thermo(name,diags):
      coeffs = [[],[]]
      with open('data/thermo.dat','r') as fin:
            lines_after_header = fin.readlines()[10:]
      for line in lines_after_header:
            lineno = lines_after_header.index(line)
            l = line.split()
            if not l:
                  continue
            elif l[0] == name:
                  if diags:
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

def calc_escape_regime(epsilon,mp,rp,T,mu,VMR,f_c_p):
      photon_energy_mass = epsilon*rp*hnu/(4*G*m_H)
      c_p_mean = sum([f_c_p(envelope_species[i],T)*VMR[i] for i in range(len(VMR))]) #J/K/mol
      gamma = c_p_mean/(c_p_mean - Rconst) #cv = cp - R; gamma = cp/cv; dimensionless
      cs2 = gamma*(Rconst/(mu/1000.))*T #remember, mu is in amu 
      alpha_B = 2.6E-19*(T/1E4)**-0.7 #cm3/s to m3/s conversion included; from Owen & Alvarez (2016)
      H_min = min(rp/3,cs2*(rp)**2./(2*G*mp))
      r_sound = G*mp/(2*cs2) #m2/s2
      if rp > r_sound:
            exit(f"R_sound = {r_sound:10.3e} < R_p = {rp:10.3e}")
      W_0_term = lambertw((-1*(rp/r_sound)**-4)*exp(3 - 4*r_sound/rp)) #dimensionless
      #NB: the lack of parentheses in eq. 16, 19, and 20 compared to Cranmer (2004)!!!
      J_cofac_inv = -4*cs2/(alpha_B*H_min) #seconds
      J_0_RR_photon = W_0_term*J_cofac_inv #photons/s
      J_0_RR_energy = W_0_term*J_cofac_inv/(4*photon_energy_mass/mp) #photon/s
      return photon_energy_mass,J_0_RR_photon,J_0_RR_energy

def calculate_water_photolysis():
      Xrayrate = 0.; EUVrate = 0.; NUVrate = 0.
      #X-ray (5-100 A), EUV (100-920 A), NUV (920-1940 A)
      #xray = [i for i in range(5,100)]; euv = [i for i in range(100,920)]; nuv = [i for i in range(920,1940)]
      water_w_mks = []; water_c_mks = []; repeats = []
      spectra = fits.getdata('data/photo/GJ581_adapt-const-res-sed.fits',1)
      spec_w = spectra['WAVELENGTH'] #Angstroms
      spec_f = spectra['FLUX'] #erg/cm2/s/Angstrom
      d_GJ581 = 6.3*206265 #GJ581 is 6.3 pc away; convert parsecs to au
      scaling = d_GJ581**2/a_planet**2 #au/au
      spec_w_mks = [w*1E-10 for w in spec_w] #convert to m
      spec_f_mks = [f*scaling*1E-3 for f in spec_f] #convert to J/m2/s/A = 1E-7*1E4
      #convert to photons/cm2/s
      spec_f_photon = [f*w/(h*c) for f,w in zip(spec_f_mks,spec_w_mks)] #photons/m2/s/A
      with open('data/photo/h2o.txt','r') as water: #data is in microns and cm2/molecule
            wlines = water.readlines()[5:]
            for line in wlines:
                  l = line.split()
                  #convert from A and cm2/molecule to m and m2/molecule while reading in
                  w = 1.E-6*float(l[0]); water_c_mks.append(1E-4*float(l[1]))
                  if w in water_w_mks:
                        w = water_w_mks[-1]*1.0000001
                  water_w_mks.append(w)
      f_h2o = interp1d(water_w_mks,water_c_mks,kind='cubic')
      f_spec = interp1d(spec_w_mks,spec_f_photon,kind='cubic')
      #water predominantly ionizes shortward of ~90 nm, rather than dissociating; index=0/10/830 at 90/100/920A, 
      convolve = [f_h2o(l*1E-10)*f_spec(l*1E-10) for l in range(90,1939)] #=m2/molecule * (photon/m2/s/A) * 1 A = /s
      Xrayrate = sum(convolve[:9])
      EUVrate = sum(convolve[10:829])
      NUVrate = sum(convolve[830:])
      return Xrayrate,EUVrate,NUVrate

def new_func():
      return
