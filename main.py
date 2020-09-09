#This contains the code to run the escape code, given the user-
#     specified values in planet.py. Note that the casual user 
#     should think before modifying the order of operations.

#Standard library imports
import matplotlib.pyplot as plt
from math import pi,log10,log,exp
from sys import exit
from numpy import linspace,logspace,argmax
from scipy.interpolate import interp1d,interp2d

#Custom functions and constants
from modules import read_baraffe,read_hidalgo
from modules import read_thermo
from modules import generic_diffusion,crossover_mass,planet_radius
from modules import find_nearest,find_above,find_below,f_lum_CW18
from modules import calculate_water_photolysis,calc_escape_regime
from atmospheric_mass_estimate import atm_est
from constants import m_Earth,r_Earth,flux_Earth,m_ocean
from constants import l_Sun,m_Sun,r_Sun,xuv_Sun
from constants import sigma,m_H,G,kb,Rconst,N_A
from constants import au2m,m2cm,yr2s,ergcm2s2Wm2
from constants import euv_a,euv_b
from planet import m_planet,r_planet,a_planet,core_frac,core_den,albedo
from planet import envelope_frac,envelope_comp,envelope_compm,envelope_species,envelope_compH
from planet import m_star,age_star,l_star,d_star,J_mag_star
from planet import stellar_tracks,norm_lum,do_euv_sat,do_emp_sat,do_emp_scale
from planet import xuv_threshold,p_photo,p_xuvcofac,mode,efficiencies,variable_efficiency,hnu

plots = True
estimate_atmosphere = False
calc_water_photo = False
diags = True
save_plots = True; plotdir = "saved_plots/"

#Read in initial luminosity data, establish finer age grid and interpolate luminosity
if stellar_tracks =='Baraffe': #Baraffe et al. (2015; A&A)
      a1,a2,l1,l2,m1,m2 = read_baraffe(m_star)
elif stellar_tracks =='Hidalgo': #Hidalgo et al. (2018; ApJ)
      a1,a2,l1,l2,m1,m2 = read_hidalgo(m_star)
if diags:
      print(f"Star is located between {m1:5.3f} and {m2:5.3f} solar masses; using {stellar_tracks} et al. luminosity evolution")

#set minimum and maximum ages for interpolation; Hidalgo et al. data is somewhat coarse, and can
#     generate sinusoidal noise at very early times if not caught and set aside. Terrestrial planets
#     like the Earth take roughly 0.03-0.1 Gyr to form; giant planets are ~0.001-0.01 Gyr (Helled et al., 2013)
min_age = max([min(a1),min(a2),0.01]); max_age = min([max(a1),max(a2),14.5])
age = [a for a in logspace(log10(min_age),log10(max_age),num=3000,endpoint=True)]
flum1 = interp1d(a1,l1,kind='cubic'); flum2 = interp1d(a2,l2,kind='cubic') 
lum1inter = flum1(age); lum2inter = flum2(age)
f1 = (m2 - m_star)/(m2 - m1); f2 = (m_star - m1)/(m2 - m1)
luminosity = [10**(f1*log10(l1)+f2*log10(l2)) for l1,l2 in zip(lum1inter,lum2inter)]
match = f"{f1*m1 + f2*m2:4.2f}*"
if norm_lum:
      Ltest = luminosity[find_nearest(age,age_star)]
      luminosity = [l*l_star/Ltest for l in luminosity]
      print(f"L_interp({age_star:4.2f} Gyr) = {Ltest:6.4f} vs. L_obs = {l_star:6.4f} L/L_sun ({100*abs(l_star-Ltest)/l_star:6.3f}% difference)")

#Identify the 'start' of the main sequence, which is roughly the minimum in the luminosity
#     This does not work for very low-mass objects, so I've set it to exit gracefully.
ms_start = -1
ms_start = age[luminosity.index(min(luminosity))]
if ms_start == max(age):
      exit("Main sequence start was not found - check data, make sure m_star>~0.08M_Sun")
current_age_ind = find_nearest(age,age_star)
current_age = age[current_age_ind]
if diags:
      print(f"Main sequence starting at {round(ms_start,4)} Gyr")
      print(f"Reported age of system: {round(age_star,3)} Gyr and {stellar_tracks} nearest age = {round(age[current_age_ind],3)} Gyr")
      lum_guess = f_lum_CW18(m_star)
      print(f"Cuntz & Wang (2018) L = {lum_guess:10.2e} vs. L(current_age) = {luminosity[current_age_ind]:10.2e}"
          +f"; difference of  ~{100*(luminosity[current_age_ind] - lum_guess)/luminosity[current_age_ind]:6.2f}%"
          +f"\n L(MS start) = {luminosity[age.index(ms_start)]:10.2e}, and L-L_CW difference is ~{100*(luminosity[age.index(ms_start)] - lum_guess)/luminosity[age.index(ms_start)]:6.2f}%")

#calculate the flux at the planet; flux and flux_Earth have units of ergs/cm2/s
flux = [l*l_Sun/(4*pi*(a_planet*au2m*m2cm)**2.) for l in luminosity]
flux_scaled = [f/flux_Earth for f in flux]
if diags:
      print(f"Present-day planet instellation = {round(flux_scaled[age.index(current_age)]*100,4)}% Earth's insolation")
      print(f"    Instellation @ 100 Myr = {round(flux_scaled[find_nearest(age,0.1)]*100,4)}% F_Earth")

#find saturation timescale
tau_sat = 8. #default of log10(0.1) Gyr

#tau_i calculation from Sanz-Forcada et al. (2011; A&A) for L_X = sum(5 A < lambda < 920 A)
tau = lambda index : 2.03E20*(luminosity[index]*l_Sun)**(-0.65)
tau_sat = age[min(range(len(age)), key = lambda i: abs(age[i]-tau(i)))]
if do_emp_sat:
      temp = [(1.89E28*a**-1.55)/(4*pi*(a_planet*au2m*m2cm)**2.) < 6.3E-4*f for a,f in zip(age,flux)]
      tau_sat = age[argmax(temp)]
if diags:
      if do_emp_sat and not do_emp_scale:
            print(f"Sanz-Forcada et al. (2011) saturation time is {tau_sat:5.3} Gyr.")
            print("        ---> Using empirically-fit 'saturation' time to match saturated/subsaturated fluxes from eq. 5")

#X-ray (5-100 A) + EUV (100-920 A) = XUV (5-920 A)
xray_flux = []; euv_flux = []; xuv_flux = []
xray_scale = []; euv_scale = []; saturation = []

if do_emp_scale:
#Peacock et al. (2020) reports F_x/F_J and F_euv/F_J, which we scale with F_J_bol:
#     F_J/F_bol = 2.5119^(m_bol - m_J), where m_bol = 5log(d_star/10 pc) - 2.5log(L/L_sun) 
#                                         and m_J is stated in TOI-1266 discovery paper
      F_J_bol = 2.5119**(5*log10(d_star/10)-2.5*log10(l_star) - J_mag_star)
      f_euv = lambda age_gyr : 2.38*(age_gyr*1e3)**-1
      f_xray = lambda age_gyr : 10.78*(age_gyr*1e3)**-1.36
      tau_euv = age[min(range(len(age)), key = lambda i: f_euv(age[i]) > 8.6E-3)]
      tau_xray = age[min(range(len(age)), key = lambda i: f_xray(age[i]) > 5.3E-3)]
      tau_sat = max(tau_euv,tau_xray)
      if diags:
            print("Using EUV/X-ray empirical fit for M dwarfs from Peacock et al. (2020)")
            print(f"    Saturation time is set by larger of EUV & X-ray saturation at {tau_sat:5.3f} Gyr.")
            print(f"    ----> F_J/F_bol = {F_J_bol:8.2e}")

for i in range(len(age)):
      if do_emp_scale: #empirical scaling from Peacock et al. (2020) - HAZMAT VI
            euv_lum = luminosity[i]*l_Sun*min(8.6E-3,f_euv(age[i]))*F_J_bol
            xray_lum = luminosity[i]*l_Sun*min(5.3E-3,f_xray(age[i]))*F_J_bol
      else:
            euv_lum = 10**(29.12+euv_a - (1.24+euv_b)*log10(age[i]))
            if do_euv_sat:
                  euv_lum = min(euv_lum, 8.6E-3*luminosity[i]*l_Sun)
            if age[i]<tau_sat:                  #saturated regime
                  xray_lum = 6.3E-4*luminosity[i]*l_Sun
            else:                               #undersaturated regime
                  xray_lum = 1.89E28*age[i]**-1.55
      xray_flux.append(xray_lum/(4*pi*(a_planet*au2m*m2cm)**2.))
      euv_flux.append(euv_lum/(4*pi*(a_planet*au2m*m2cm)**2.))
      xuv_flux.append(euv_flux[i]+xray_flux[i])
      #scale these values to make plotting a little easier
      saturation.append(xuv_flux[i]/flux[i])
      euv_scale.append(euv_flux[i]/flux[i])
      xray_scale.append(xray_flux[i]/flux[i])
xuv_flux_Earth = [x/xuv_Sun for x in xuv_flux]

#Handy lambda functions for dependent variables
f_roche = lambda rp,mp,ms,a : a*au2m*((mp*m_Earth)/(3*ms*m_Sun))**(1/3)/(rp*r_Earth)      #in planet radii
f_ktide = lambda ro,rp : (1 - (3/(2*(ro/rp))) + (1/(2*(ro/rp)**3)))                       #dimensionless
f_T_eq  = lambda f,a : ((f*(1-a))/(4*sigma))**0.25
f_grav  = lambda mp,rp : (G*mp*m_Earth)/(rp*r_Earth)**2.
f_mu    = lambda xi,mi : 1/sum(xi[i]/mi[i] for i in range(len(xi))) #requires mass fraction as m; in amu
f_vmr   = lambda mi,mu,mass : [mi[i]*mu/mass[i] for i in range(len(mi))] #mass fraction m
f_mult  = lambda xi,ei : sum(xi[i]*ei[i] for i in range(len(xi))) #generic vector multiplication
f_H_below=lambda t,mu,g : kb*t/(mu*m_H*g)
f_rho_P = lambda pressure,temperature,mu : [mu*m_H*p/(kb*t) for p,t in zip(pressure,temperature)]
f_xuvcr = lambda mp,rp,eff: 180*mp**2*rp**-3*(eff/0.3)**-1 #Critical XUV flux for O escape from Luger & Barnes (2018)

#Prepare for call to escape regime by reading in thermo data for each species and defining the equation for them
coeffs = {}
for species in envelope_species:
      coeffs[species] = read_thermo(species,diags)
#c_p = [J/K/mol]
f_c_p = lambda n,T : ((T<1000)*(sum([coeffs[n][0][i+2]*T**i for i in range(-2,5)])) + (T>=1000 and T<6000)*(sum([coeffs[n][1][i+2]*T**i for i in range(-2,5)])))

#define initial variable states
c_f_t = [core_frac]; c_d_t = [core_den]
e_f_t = [envelope_frac]; e_comp_t = [envelope_comp]; num_comps = len(envelope_comp)
m_p_t = [m_planet]
r_p_t = [planet_radius(m_p_t[-1],0,100)] #[r_planet] Can't use observed radius in combination with M/R relationship below
grav_t = [f_grav(m_planet,r_planet)] #Gravity [m/s2] at the 20 mbar level assumed to be for m_planet @ r_planet
roche_t = [f_roche(r_p_t[-1],m_p_t[-1],m_star,a_planet)]
ktide_t = [f_ktide(roche_t[-1],r_p_t[-1])]
mu_t = [f_mu(e_comp_t[-1],envelope_compm)]
vmr_t = [f_vmr(e_comp_t[-1],mu_t[-1],envelope_compm)]
Htot = sum([v*H for v,H in zip(vmr_t[-1],envelope_compH)]) #assuming that H in other species is available for escape
epsilon_t = [max(efficiencies[0]*Htot,f_mult(vmr_t[-1],efficiencies))]#f_mult(vmr_t[-1],efficiencies)] 
escape_flux_t = []; mescape_flux_t = []
crossover_mass_t = []; mass_change = 0.
new_mass = [0 for i in range(len(envelope_species))]
F_RR2photon = []; F_RR2energy = []; mp_crit = []
RR_loss = 0.; energy_loss = 0.; photon_loss = 0.; diffusion_loss = 0.
regime_t = []
critical_xuv_time = -1.

if calc_water_photo:
      #calculate how much water can be photolyzed through the star's age, assuming there's only H2O
      h2o_x,h2o_euv,h2o_nuv = calculate_water_photolysis()
      TA_temp = f_T_eq(flux[current_age_ind]*ergcm2s2Wm2,albedo)/(4*pi*(r_planet*r_Earth)**2) #assume layer is 1 m deep
      #/s * molecules/m3 * m3 * kg/moles * moles/molecules = kg/s
      h2o_x_mass = h2o_x * (p_photo/(kb*TA_temp)) * 0.018 / N_A 
      h2o_euv_mass = h2o_euv * (p_photo/(kb*TA_temp)) * 0.018 / N_A
      h2o_nuv_mass = h2o_nuv * (p_photo/(kb*TA_temp)) * 0.018 / N_A
      if diags:
            print(f"GJ581 present-day spectra can photolyze {h2o_x_mass+h2o_euv_mass+h2o_nuv_mass:8.3e} kg/s")

#Establish the bounds for the modes
if mode == 'forward':
      start = [0]
      end = [len(age)]
      step = [1]
      print("Solving for initial planet mass from starting mass,radius, and composition in forward mode.")
elif mode == 'reverse':
      start = [current_age_ind,current_age_ind+1]
      end = [-1,len(age)]
      step = [-1,1]
      print("Solving for initial planet mass from present mass, radius, and composition in reverse mode.")
      if e_comp_t[0] == 0.:
            exit("Reverse mode logic is incompatible with no present-day H2 at this time. Try again tomorrow.")
else:
      exit("Unrecognized mass/radius/composition mode.")

#Have to define this in here - T_eq depends on when we start.
T_eq_t = [f_T_eq(flux[start[0]]*ergcm2s2Wm2,albedo)]
      
if diags:
      print(f"Starting from {round(age[start[0]],3)} Gyr, L/L_sol = {luminosity[start[0]]:6.4f}, XUV/L_bol = {saturation[start[0]]:8.2e}")
      print('Starting values for some variables:')
      print('Roche = ',roche_t,' Ktide = ',ktide_t,' mu = ',mu_t,' T_eq = ',T_eq_t,' epsilon = ',epsilon_t)
      print('core_frac = ',c_f_t,' envelope_frac = ',e_f_t)
      print(f' envelope_comp: {envelope_species} = {envelope_comp}')

#Estimating the atmospheric mass relies on the equilibrium temperature
if estimate_atmosphere:
      T_skin = T_eq_t[-1]/2**0.25
      fe_core_radius = planet_radius(c_f_t[-1]*m_p_t[-1],100,0)
      h2o_core_radius = planet_radius(c_f_t[-1]*m_p_t[-1],0,100)
      atm_est(grav_t[-1],T_skin,e_f_t[-1]*m_p_t[-1]*m_Earth,fe_core_radius*r_Earth,vmr_t[-1],mu_t[-1],diags,plotdir)
      if diags:
            print(f"Fe core = {fe_core_radius:6.2f} Earth radii; water core = {h2o_core_radius:6.2f} Earth radii")
            
      
#Now solve for the planet mass,radius, and composition through time
#Variables to update: core_frac/den, envelope_frac/comp/compm, mu, m_planet, r_planet, roche, ktide,
#                       escape_flux,crossover_mass,T_eq, 
#Order of operations: 1) calculate escape regime; 2) calculate escape flux; 3) calculate crossover_mass
#                       4) update mass into relevant bins
for j in range(len(start)):
      #now we do the time integration
      for i in range(start[j],end[j],step[j]):
            #define common variables for both regimes
            if i > 0:
                  dts = (age[i]-age[i-1])*1E9*yr2s
            else:
                  dts = (age[1]-age[0])*1E9*yr2s
            H_below = f_H_below(T_eq_t[-1],mu_t[-1],grav_t[-1])         #m
            p_xuvbase = p_xuvcofac*grav_t[-1]                           #Pa
            r_p_m = r_p_t[-1]*r_Earth                                   #m
            rbase = r_p_m + H_below*log(p_photo/p_xuvbase)              #m
            m_p_kg = m_p_t[-1]*m_Earth                                  #kg
            #TODO Approximate the flux-dependent mass loss efficiency of Bolmont et al. (2017), Fig. 2
            #Koskinen et al. (2014) efficiencies for ultra-HJs are 0.085@0.2 au; 0.22@0.1 au; 0.44@<0.1 au
            #if scale_efficiency: 
            #
            #Start by determining which regime the escape occupies by calculating J_0 from eqns. 18-20 in Owen & Alvarez (2016)
            #Returning the threshold planetary mass for the energy- to photon-limited regimes (Mp > actual mass means photon-limited),
            #     and the flux (in photon/m2/s) needed to hop between the different regimes.
            if xuv_flux[i] > 180.: #Kelvin - motivated by Koskinen et al. (2007) and Murray-Clay et al. (2009)
                  T_exo = 2E4 #Alternatively, 2E4, but <15% change in value
                  #gamma_exo = 1.53 #Determined from McBride et al. for H, H2, and O = 1.66 @ 10k K -> 1.53 @ 20k K
            else:
                  T_exo = 2E3 
                  #c_p_exo = sum([f_c_p(envelope_species[i],T_exo)*vmr_t[-1][i] for i in range(len(vmr_t[-1]))])*Rconst #J/K/mol
                  #gamma_exo = c_p_exo/(c_p_exo - Rconst)
            Mp_p2e, flux1, flux2 = calc_escape_regime(epsilon_t[-1],m_p_kg,r_p_m,T_exo,mu_t[-1])
            #convert both fluxes from photon/m2/s to erg/cm2/s
            flux1 = flux1*hnu/ergcm2s2Wm2
            flux2 = flux2*hnu/ergcm2s2Wm2
            mp_crit.append(Mp_p2e/m_Earth)
            F_RR2photon.append(flux1.real)
            F_RR2energy.append(flux2.real)
            if critical_xuv_time < 0. and xuv_flux[i]<f_xuvcr(m_p_t[-1],r_p_t[-1],epsilon_t[-1]):
                  critical_xuv_time = age[i]
                  if diags:
                        print(f"Critical XUV flux (O loss) = {f_xuvcr(m_p_t[-1],r_p_t[-1],epsilon_t[-1]):8.2e} erg/cm2/s")
                        print(f'\nPlanet drops below L&B F_XUV_crit at {critical_xuv_time:6.4f} Gyr.')
                        print(f"XUV flux falls below 40x Earth's (~{xuv_flux[find_nearest(xuv_flux_Earth,40)]:6.2f} erg/cm2/s) at {age[find_nearest(xuv_flux_Earth,40)]:6.4f} Gyr.\n")
            if diags and i == current_age_ind:
            #      print(epsilon_t[-1],m_p_kg,r_p_m,T_eq_t[-1],mu_t[-1],vmr_t[-1])
                  print(f'{age[i]:5.3f}-Gyr critical mass = {mp_crit[-1]:8.3f} Earth masses\n'+
                        f'RR/photon flux threshold = {F_RR2photon[-1]:10.2e} erg/cm2/s\n'+
                        f'RR/energy flux threshold = {F_RR2energy[-1]:10.2e} erg/cm2/s\n')
            #include diffusion limit
            ref_H_flux = epsilon_t[-1]*xuv_flux[i]*ergcm2s2Wm2*r_p_m/(4*G*m_p_kg*ktide_t[-1]*m_H) #H/m2/s
            Htot = sum([v*H for v,H in zip(vmr_t[-1],envelope_compH)]) #assuming that H in other species is available for escape
            cm,diff_limit = crossover_mass(T_eq_t[-1],ref_H_flux,grav_t[-1],Htot,mu_t[-1])
            crossover_mass_t.append(cm)
            #Select escape regime based on XUV threshold fluxes
            if m_p_t[-1] < mp_crit[-1] and xuv_flux[i] < F_RR2photon[-1]:
                  #Photon-limited escape rate, per Owen & Alvarez (2016)
                  mflux = pi*r_p_m**2*m_H*xuv_flux[i]*ergcm2s2Wm2/hnu #convert from W/m2 to photons/m2/s; kg/s
                  regime_t.append(0)
                  photon_loss = photon_loss + mflux*dts/m_Earth
            elif xuv_flux[i] > F_RR2energy[-1]:
                  #calculate escape rate with radiation-recombination limit from Murray-Clay et al. (2009, ApJ) as 
                  #     normalized by Luger et al. (2015, Astrobiology)
                  mflux = 7.11E4*(xuv_flux[i])**0.5*r_p_t[-1]**1.5 #kg/s
                  regime_t.append(2)
                  RR_loss = RR_loss + mflux*dts/m_Earth
            elif Htot > 0.01: #If [H] is not prevalent, we enter the diffusion limited regime
                  #calculate energy-limited escape rate following Luger & Barnes (2015)
                  mflux = epsilon_t[-1]*pi*xuv_flux[i]*ergcm2s2Wm2*r_p_m*rbase**2/(G*m_p_kg*ktide_t[-1]) #kg/s
                  regime_t.append(1)
                  energy_loss = energy_loss + mflux*dts/m_Earth
            else:
                  mflux = diff_limit*4*pi*r_p_m**2/m_H
                  regime_t.append(4)
                  diffusion_loss = diffusion_loss + mflux*dts/m_Earth
            eflux = mflux/(m_H*4*pi*r_p_m**2.) #H atoms/m2/s; assume the escape flux is H to start, then used modified drag-off below
            mescape_flux_t.append(mflux)
            escape_flux_t.append(eflux)
            #Convert escape flux into planet mass and inventory modifications
            #Include only species that are below the crossover mass and have non-zero inventories
            #TODO what about drag-off in reverse mode? It should only see if a gas is already present, otherwise no drag-off or re-addition?
            escaping_species = [envelope_compm.index(m) for m in envelope_compm[1:] if crossover_mass_t[-1]>m and e_comp_t[-1][envelope_compm.index(m)]>0.] 
            new_mass = [e_f_t[-1]*m_p_kg*e_comp_t[-1][s] for s in range(len(envelope_species))]
            old_mass = new_mass #for differencing later, once all the escape logic is carried out
            delta_mass = 0.
            if any(escaping_species):
                  escaping_MMRs = [max(1.E-30,e_comp_t[-1][s]) for s in escaping_species]
                  escaping_MMRs.insert(0,max(1.E-30,e_comp_t[-1][0]))
                  escaping_VMRs = [max(1.E-50,vmr_t[-1][s]) for s in escaping_species]
                  escaping_VMRs.insert(0,Htot)
                  escaping_mus = [envelope_compm[s] for s in escaping_species]
                  escaping_mus.insert(0,1.)
                  mu_heavy = f_mu(escaping_MMRs[1:],escaping_mus[1:])*sum(escaping_MMRs[1:])
                  vmr_heavy = sum([vmr_t[-1][s] for s in escaping_species])
                  mu_escaping = f_mu(escaping_MMRs,escaping_mus)*sum(escaping_MMRs)
                  b_H = generic_diffusion(1.,mu_t[-1],T_eq_t[-1])
                  phi_prime = kb*T_eq_t[-1]*eflux/(b_H*grav_t[-1]*mu_escaping*m_H)
                  m_c_prime = 1 + ((mu_heavy - 1)*mu_heavy*vmr_heavy)/mu_escaping + phi_prime
                  crossover_mass_t[-1] = m_c_prime
                  escaping_species.insert(0,0) #need to include H now in escape
                  F_H_prime = eflux*Htot*(m_c_prime - 1)/sum([escaping_mus[k]*escaping_VMRs[k]*(m_c_prime - escaping_mus[k]) for k in range(len(escaping_species))])
                  for s in escaping_species:
                        k = escaping_species.index(s)
                        del_m_comp = 4*pi*r_p_m**2*step[j]*dts*envelope_compm[s]*m_H*F_H_prime*escaping_VMRs[k]*(m_c_prime - escaping_mus[k])/(Htot*(m_c_prime - 1))
                        new_mass[s] -= del_m_comp
                        delta_mass += del_m_comp
            #TODO new mass addition based on composition, crossover mass !!!!!!!!
            else:
                  new_mass[0] -= mflux*step[j]*dts #kg
                  delta_mass += mflux*step[j]*dts #kg
            #If a mass reservoir that's supposed to be dragged off goes negative, that extra mass is lost as H
            new_mass[0] = new_mass[0] - abs(sum([m for m in new_mass if m < 0]))
            if new_mass[0]<=0 and new_mass[envelope_species.index('H2O')] > 0: #H2 exhausted; accumulate O from water
                  ind_H2O = envelope_species.index('H2O')
                  delta_O = abs(new_mass[0])*envelope_compm[ind_H2O]/envelope_compm[envelope_species.index('H2')]
                  new_mass[ind_H2O] -= delta_O
                  new_mass[envelope_species.index('O')] += delta_O*(16/18)
                  #TODO finish this implementation - small impact, since O accumulates and then escapes in next step?
                  #new_F_H = abs(new_mass[0])/(m_H*dts*4*pi*r_p_m**2) #convert H lost from water into H flux to see if it drags off O
                  #if crossover_mass(T_eq_t[-1],new_F_H,grav_t[-1],Htot,mu_t[-1]) > 16:
                  #      F_H_prime = eflux*(2/3)*vmr_t[-1][ind_H2O]*(m_c_prime - 1)/(envelope_compm[ind_H2O]*vmr_t[-1][ind_H2O]*(m_c_prime - envelope_compm[ind_H2O]))
                  new_mass[0] = 0.
            new_mass = [max(0,m) for m in new_mass] #at the end, if H2 is gone and water is depleted, those reservoirs should be zero
            #delta_mass = -1*sum([old-new for old,new in zip(old_mass,new_mass)])
            new_e_frac = [m/max(1.E-90,sum(new_mass)) for m in new_mass]
            mass_change = mass_change+abs(delta_mass)/m_Earth
            #if crossover_mass_t[-1] > 18:
            #      print(age[i],regime_t[-1],epsilon_t[-1],m_p_kg,r_p_m,mu_t[-1],vmr_t[-1])
            #else:
            #      print('nah -- ',regime_t[-1],epsilon_t[-1],m_p_kg,r_p_m,mu_t[-1],vmr_t[-1])
            
            #update variables
            if i < len(age)-1:
                  #Append new values as we move backwards in time, then 'flip' variables
                  m_p_t.append(m_p_t[-1]-delta_mass/m_Earth)
                  c_f_t.append(core_frac*m_planet/m_p_t[-1])
                  e_f_t.append(1.-c_f_t[-1])
                  e_comp_t.append(new_e_frac)
                  r_p_t.append(planet_radius(m_p_t[-1],0,100)) #Assuming that the planet follows the 100% H2O M-R relation from Noack et al. (2016)
                  grav_t.append(f_grav(m_p_t[-1],r_p_t[-1]))
                  roche_t.append(f_roche(r_p_t[-1],m_p_t[-1],m_star,a_planet))
                  ktide_t.append(f_ktide(roche_t[-1],r_p_t[-1]))
                  T_eq_t.append(f_T_eq(flux[i]*ergcm2s2Wm2,albedo))
                  if any([m>0 for m in e_comp_t[-1]]):
                        mu_t.append(f_mu(e_comp_t[-1],envelope_compm))
                        vmr_t.append(f_vmr(e_comp_t[-1],mu_t[-1],envelope_compm))
                        epsilon_t.append(max(efficiencies[0]*Htot,f_mult(vmr_t[-1],efficiencies)))
                  else:
                        mu_t.append(44)
                        vmr_t.append([0 for s in envelope_species])
                        epsilon_t.append(0.01)
      
      if mode == 'reverse' and j == 0:
            #Flip variables so that they're aligned with *age* variable
            c_f_t.reverse(); c_d_t.reverse()
            e_f_t.reverse(); e_comp_t.reverse()
            r_p_t.reverse(); m_p_t.reverse()
            mu_t.reverse(); T_eq_t.reverse()
            roche_t.reverse(); ktide_t.reverse()
            mp_crit.reverse(); F_RR2photon.reverse(); F_RR2energy.reverse()
            regime_t.reverse()
            epsilon_t.reverse(); crossover_mass_t.reverse()
            mescape_flux_t.reverse()
            escape_flux_t.reverse()
            mp_crit.reverse(); F_RR2photon.reverse(); F_RR2energy.reverse()

if diags:
      m_oxygen_gained = e_f_t[current_age_ind]*m_p_t[current_age_ind]*e_comp_t[current_age_ind][envelope_species.index('O')]
      print(f"Photon-limited mass loss = {photon_loss:10.4e}; RR-limited loss = {RR_loss:10.4e}; energy-limited loss = {energy_loss:10.4e}; diffusion-limited loss = {diffusion_loss:10.4e} [m_Earth]")
      print(f"{m_p_t[0]:6.4f}-->{m_p_t[-1]:6.4f} m_Earth; total mass change = {mass_change:6.4f} m_Earth")
      print(f"{r_p_t[0]:6.4f}-->{r_p_t[-1]:6.4f} r_Earth over this time, or {100*(r_p_t[-1]-r_p_t[1])/r_p_t[0]:6.4f}%")
      print(f"Present-day escape rate = {mescape_flux_t[age.index(current_age)]:8.2e} kg/s")
      print(f"                        = {escape_flux_t[age.index(current_age)]:8.2e} molec/m2/s")
      print(f"                        = {escape_flux_t[age.index(current_age)]/1e4:8.2e} molec/cm2/s")
      print(f"\n  M_O2 gained = {m_oxygen_gained:7.5f} Earth masses ({m_oxygen_gained*m_Earth/(16*m_ocean/18):7.5f} Earth oceans equiv.) at current age")
      print(f"Present-day m_p = {m_p_t[current_age_ind]:6.4f} m_Earth; r_p = {r_p_t[current_age_ind]:6.4f} r_Earth")
      if mode == 'reverse':
            print(f"Initial envelope fraction = {e_f_t[0]:10.5e} and composition = {e_comp_t[0]}; initial R = {r_p_t[0]:8.6f} and M = {m_p_t[0]:8.6f}")
      else:
            print(f"Present-day envelope fraction = {e_f_t[current_age_ind]:10.5e} and composition = {e_comp_t[current_age_ind]}; initial R = {r_p_t[current_age_ind]:8.6f} and M = {m_p_t[current_age_ind]:8.6f}")
      print(f"Final envelope fraction = {e_f_t[-1]:10.5e} and composition = {e_comp_t[-1]}; R = {r_p_t[-1]:8.6f} and M = {m_p_t[-1]:8.6f}")
if plots:
      fig1,ax1 = plt.subplots()
      if diags:
            ax1.plot(a1,l1,'r--') #lower stellar mass for luminosity interpolation
            ax1.plot(a2,l2,'r--') #upper stellar mass for luminosity interpolation
#            ax1.plot([1.1*max(age),1.1*max(age)],[0.98*f_lum_CW18(m_star),1.02*f_lum_CW18(m_star)],color='gray')
      ax1.plot(age,luminosity,'k')
      #include saturation timescale
      ax1.plot(tau_sat,luminosity[age.index(tau_sat)],'bo')
      ax1.text(tau_sat,1.1*luminosity[age.index(tau_sat)],r'${\tau}_{sat}$',color='b',ha='right')
      #include the notional start of the main sequence lifetime for star
      ax1.plot(ms_start,luminosity[age.index(ms_start)],'ro',markerfacecolor='None')
      ax1.text(ms_start,1.1*luminosity[age.index(ms_start)],'MS',color='r',ha='left')
      #include the notional age of the system
      ax1.plot(current_age,luminosity[current_age_ind],'x',color='gray')
      ax1.text(current_age,1.1*luminosity[current_age_ind],'Today',color='gray',ha='center')
      ax1.set_xscale('log'); ax1.set_yscale('log')
      ax1.set_xlabel("Age [Gyr]"); ax1.set_ylabel("Luminosity [L/L$_{\odot}$]")
      ax1.set_xlim(0.9*min(age),1.2*max(age))
      ax1.set_ylim(0.9*min(l1),1.2*max(luminosity))
      ax1.set_yticks([0.01,0.1])
      ax1.set_yticklabels(['0.01','0.1'])
#      ax1.set_title(stellar_tracks+" et al. luminosity evolution for "+str(match)+" M$_{\odot}$ star\n"
#            +"Luminosities for bookend masses(r--) and interpolation(k)")
      fig1.tight_layout()
      plt.show()
      if save_plots:
            fig1.savefig(plotdir+'luminosity_MS_tau_i_vs_time.png',dpi=200)

if plots:
      fig2 = plt.figure()
      ax2 = plt.subplot(211)
      ax2.plot(age,m_p_t,'k',label=r'Planet Mass')
      ax2.set_xscale('log')
      ax2.text(age[0],m_p_t[0],'Planet',c='k',ha='left',va='top')
      ax2.plot(age,[i*j for i,j in zip(m_p_t,e_f_t)],'b',label='Envelope Mass')
      ax2.text(age[0],m_p_t[0]*e_f_t[0],'Envelope',c='b',ha='left',va='bottom')
      ax2.plot(age,[i*j for i,j in zip(m_p_t,c_f_t)],'r',label='Core Mass')
      ax2.text(age[0],m_p_t[0]*c_f_t[0],'Core',c='r',ha='left',va='top')
      ymin, ymax = ax2.get_ylim()
      ax2.plot([current_age,current_age],[0.,ymax],':',color='gray')
      ax2.set_xlim(0.9*min(age),1.2*max(age))
      ax2.set_ylim([0.,ymax])
      ax2.set_ylabel(r"Mass [$M_{\oplus}$]")
      plt.setp(ax2.get_xticklabels(),visible=False)
#      ax2.text(current_age,0.85*m_p_t[current_age_ind],'Today',color='gray',ha='center')
      ax2b = plt.subplot(212,sharex=ax2) #ax2.twinx()
      for i in range(len(envelope_species)):
            e_comp_t_spec = [e_comp_t[j][i] for j in range(len(age))] #vmr_t[j][i] for j in range(len(age))]
            #ax2b.plot(age,e_comp_t_spec,label=envelope_species[i])
            ax2b.plot(age,[a*b*c for a,b,c in zip(m_p_t,e_f_t,e_comp_t_spec)],'--',label=envelope_species[i])
      #ax2b.set_ylabel(r'Envelope composition [VMR]')
      ax2b.set_xscale('log')
      ax2b.set_ylabel(r'Envelope composition [$M_{\oplus}$]')
      ax2b.set_xlabel('Age [Gyr]')
      ymin, ymax = ax2b.get_ylim()
      ax2b.plot([current_age,current_age],[ymin,ymax],':',color='gray')
      plt.legend(loc='center left')
      fig2.tight_layout()
      plt.show()
      if save_plots:
            fig2.savefig(plotdir+'Planet_mass_vs_time.png',dpi=200)

if plots:
      fig3,ax3 = plt.subplots()
      ax3.plot(age,saturation,'k--',label='XUV')
      ax3.plot(age,xray_scale,'r',label='X-ray')
      ax3.plot(age,euv_scale,'b',label='EUV')
      ax3.set_xscale('log'); ax3.set_yscale('log')
      ax3.set_xlabel("Age [Gyr]"); ax3.set_ylabel(r"L$_{\lambda}$/L$_{bol}$")
#      ax3.set_title(stellar_tracks+" et al.-based X-ray, EUV, and XUV evolution for "+str(match)+" M$_{\odot}$ star")
      ax3.set_xlim(0.9*min(age),1.2*max(age))
      plt.legend(loc='best')
      fig3.tight_layout()
      plt.show()
      if save_plots:
            fig3.savefig(plotdir+'UV_scaling_vs_time.png',dpi=200)
if plots:
      fig4,ax4 = plt.subplots()
      ax4.plot(age,T_eq_t,'k')
      ax4.set_xscale('log')
      ax4.set_xlabel('Age [Gyr]'); ax4.set_ylabel('Equilibrium temperature [K]')
      fig4.tight_layout()
      plt.show()
      if save_plots:
            fig4.savefig(plotdir+'equilibrium_temp_vs_time.png',dpi=200)

if plots:
      fig5,ax5 = plt.subplots()
      ax5.plot(age,xuv_flux,'k--',label='XUV')
      ax5.plot(age,xray_flux,'r',label='X-ray')
      ax5.plot(age,euv_flux,'b',label='EUV')
      plt.legend(loc='best')
      ax5.set_xscale('log'); ax5.set_yscale('log')
      ymin, ymax = ax5.get_ylim()
      ax5.set_xlim(0.9*min(age),1.2*max(age))
      ax5b = ax5.twinx()
      ax5b.loglog(age,xuv_flux_Earth,'None')
      ax5b.set_ylabel(r"XUV flux [F/F$_{\oplus}$]")
#      ax5.plot([min(age),max(age)],[1e4,1e4],color='gray',linestyle='-.')
#      print(f"XUV flux drops below 1E4 erg/cm2/s at {age[find_nearest(xuv_flux,1E4)]} Gyr.")
      ax5.set_xlabel("Age [Gyr]"); ax5.set_ylabel(r"Flux [ergs/cm2/s]")
      ax5b.set_ylim(ymin/xuv_Sun, ymax/xuv_Sun)
#      ax5.set_title(stellar_tracks+" et al.-based X-ray, EUV, and XUV evolution "+str(match)+" M$_{\odot}$ star")
      fig5.tight_layout()
      plt.show()
      if save_plots:
            fig5.savefig(plotdir+'UV_flux_vs_time.png',dpi=200)

if plots:
      fig6,ax6 = plt.subplots()
      ax6.plot(age,crossover_mass_t,'b-')
      ax6.set_xscale("log"); ax6.set_yscale('log')
      ax6.set_xlabel("Age [Gyr]"); ax6.set_ylabel('Crossover mass [amu] (b)')
      ax6.set_xlim(0.9*min(age),1.2*max(age))
      ax6b = ax6.twinx()
      ls = [':','-','--','-.']
      regime_labels = ['Photon-limited','Energy-limited','RR-limited','Diffusion-limited']
      ages = {}; fluxes = {}
      for i in range(4):
            new_mescape = [mescape_flux_t[k] if regime_t[k] == i else 'NaN' for k in range(len(age))]
            ax6b.plot(age,new_mescape,ls=ls[i],color='k',label=regime_labels[i])
      ax6b.set_yscale('log')
      #ax6b.plot(age,regime_t,'.')
      ax6b.set_ylabel('Escape flux [kg/s] (k)')
      plt.legend(loc='best')
      if min(crossover_mass_t) < 18:
            crit_mass = find_nearest(crossover_mass_t,18.)
            ax6.plot(age[crit_mass],crossover_mass_t[crit_mass],'bo')
            ax6.text(age[crit_mass+5],1.1*crossover_mass_t[crit_mass],"Water")
      fig6.tight_layout()
      plt.show()
      if save_plots:
            fig6.savefig(plotdir+'crossover_mass_vs_time.png',dpi=200)

if plots:
      fig7,ax7 = plt.subplots(2,1,sharex=True)
      ax7[0].semilogx(age,mp_crit)
      ax7[1].loglog(age,F_RR2photon,label='RR-photon')
      ax7[1].loglog(age,F_RR2energy,label='RR-energy')
      ax7[0].set_xlabel('Age [Gyr]')
      ax7[0].set_ylabel(r'Critical mass [M$_{\oplus}$]')
      ax7[1].set_ylabel(r'Flux threshold [erg/cm$^2$/s]')
      ax7[0].set_xlim(0.9*min(age),1.2*max(age))
      ax7[1].legend()
      fig7.tight_layout()
      plt.show()
      if save_plots:
            fig7.savefig(plotdir+'flux_regime_boundaries_vs_time.png',dpi=200)
