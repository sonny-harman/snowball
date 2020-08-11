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
from modules import crossover_mass,read_baraffe,read_baraffe_grid
from constants import m_Earth,r_Earth,flux_Earth
from constants import l_Sun,m_Sun,r_Sun
from constants import sigma,m_H,G,kb,Rconst
from constants import au2m,m2cm,yr2s,ergcm2s2Wm2
from constants import euv_a,euv_b
from planet import m_planet,r_planet,a_planet,core_frac,core_den,albedo
from planet import envelope_frac,envelope_comp,envelope_compm,envelope_species,mubar
from planet import m_star,age_star,l_star,d_star,J_mag_star
from planet import interp_lum,do_euv_sat,do_emp_sat,do_emp_scale
from planet import xuv_threshold,p_photo,p_xuvcofac,mode,efficiencies
from atmospheric_mass_estimate import atm_est

plots = True
estimate_atmosphere = False#True
diags = True
finediags = False
save_plots = True; plotdir = "saved_plots/"

#Handy lambda for finding the nearest value in a list.
find_nearest = lambda vector,value : min(range(len(vector)), key = lambda i: abs(vector[i]-value))

#lambda function for finding the luminosity (in solar) for a 0.2<M<0.85 M_solar star (Cuntz & Wang, 2018)
f_lum_CW18 = lambda M : M**(-141.7*M**4. + 232.4*M**3. - 129.1*M**2. + 33.29*M + 0.215)

#Read in initial Baraffe et al. (2015; A&A), establish finer age grid and interpolate luminosity
match1,a1,l1,message1 = read_baraffe(m_star)
print("Remember, if you want to do the upper end of the uncertainty in age, you'll need to extrapolate!")
#experimental two-dimensional cubic spline for fitting both mass and age
mall,ageall,tempall,lumall,gravall,limits = read_baraffe_grid()
masses = sorted(list(set(mall)))
next_mass = masses[masses.index(match1)+1]
print(f"masses = {masses}")
f_star = interp2d(ageall,mall,lumall,kind='cubic',fill_value=-1)
#new_age = logspace(log10(limits[0]),log10(limits[1]),2000,endpoint=True)
new_age = linspace(limits[0],limits[1],2000,endpoint=True)
print(f"Limits of interpolation space = [{limits[0]:8.5f}-{limits[1]:8.5f}]")
m_vec = [m_star]*len(new_age)
flum1 = interp1d(a1,l1,kind='cubic')
age1 = linspace(min(a1),max(a1),num=2000,endpoint=True)
lum1 = flum1(age1)
age1 = [a for a in age1]
zage1 = [0] + age1
lum1.tolist()
match2,a2,l2,message2 = read_baraffe(next_mass)
flum2 = interp1d(a2,l2,kind='cubic')
age2 = linspace(min(a2),max(a2),num=2000,endpoint=True)
lum2 = flum2(age2)
age2 = [a for a in age2]
zage2 = [0] + age2
lum2.tolist()
luminosity = lum1; age = age1; zage = zage1; match = match1; message = message1

if interp_lum:
      ageinter = [a for a in linspace(max(min(a1),min(a2)),min(max(a1),max(a2)),num=2000,endpoint=True)]
      zageinter = [0] + ageinter
      lum1inter = flum1(ageinter); lum2inter = flum2(ageinter)
      f1 = (match2 - m_star)/(match2 - match1); f2 = (m_star - match1)/(match2 - match1)
      print(f1,f2)
      flux1 = [l*l_Sun/(4*pi*(a_planet*au2m*m2cm)**2) for l in lum1]
      flux2 = [l*l_Sun/(4*pi*(a_planet*au2m*m2cm)**2) for l in lum2]
      luminter = [10**(f1*log10(l1)+f2*log10(l2)) for l1,l2 in zip(lum1inter,lum2inter)]
      Ltest = luminter[find_nearest(ageinter,age_star)]
      if diags:
            print(f"L_interp({age_star:4.2f} Gyr) = {Ltest:6.4f} vs. L_obs = {l_star:6.4f} L/L_sun ({100*abs(l_star-Ltest)/l_star:6.3f}% difference)")
      luminter = [l*l_star/Ltest for l in luminter]
      fluxinter = [l*l_Sun/(4*pi*(a_planet*au2m*m2cm)**2) for l in luminter]
      luminosity = luminter; age = ageinter; zage = zageinter; match = f"{f1*match1 + f2*match2:4.2f}*"; message = f"Interpolating between {match1} and {match2} solar mass stars!!!"
figx,axx = plt.subplots()
axx.scatter(ageall,mall,c=[log10(l) for l in lumall],cmap='coolwarm')
cmap_obj = plt.cm.get_cmap('coolwarm')
rgba_new = cmap_obj([log10(l)*log10(max(lumall))/log10(max(luminosity)) for l in luminosity])
axx.scatter(age,[m_star]*len(age),c=rgba_new)
#plt.show()
#exit()

#Identify the 'start' of the main sequence, which is roughly the minimum in the luminosity
#     This does not work for very low-mass objects, or potentially for those with humps in
#     luminosity evolution. Would check twice for any other cases.
ms_start = -1
for i in range(len(age)-2):
      if luminosity[i+1]>=luminosity[i] and luminosity[i+1]<luminosity[i+2]:
            ms_start = age[i+1]
            break
if ms_start<0:
      exit("Main sequence start was not found - check data, make sure m_star>~0.08M_Sun")
current_age_ind = find_nearest(age,age_star)
current_age = age[current_age_ind]
if diags:
      print(message)
      print(f"Main sequence starting at {round(ms_start,4)} Gyr")
      print(f"Reported age of system: {round(age_star,3)} Gyr and Baraffe nearest age = {round(current_age,3)} Gyr")
      lum_guess = f_lum_CW18(m_star)
      print(f"Cuntz & Wang (2018) L = {lum_guess:10.2e} vs. L(current_age) = {luminosity[current_age_ind]:10.2e}"
          +f"; difference of  ~{100*(luminosity[current_age_ind] - lum_guess)/luminosity[current_age_ind]:6.2f}%"
          +f"\n L(MS start) = {luminosity[age.index(ms_start)]:10.2e}, and L-L_CW difference is ~{100*(luminosity[age.index(ms_start)] - lum_guess)/luminosity[age.index(ms_start)]:6.2f}%")

#calculate the flux at the planet; flux and flux_Earth have units of ergs/cm2/s
flux = [l*l_Sun/(4*pi*(a_planet*au2m*m2cm)**2.) for l in luminosity]
flux_scaled = [f/flux_Earth for f in flux]
if diags:
      print(f"Present-day planet instellation = {round(flux_scaled[age.index(current_age)]*100,4)}% Earth's insolation")

#find saturation timescale
tau_sat = 8. #default of log10(0.1) Gyr

#tau_i calculation from Sanz-Forcada et al. (2011; A&A) for L_X = sum(5 A < lambda < 920 A)
tau = lambda index : 2.03E20*(luminosity[index]*l_Sun)**(-0.65)
tau_sat = age[min(range(len(age)), key = lambda i: abs(age[i]-tau(i)))]
if do_emp_sat:
      temp = [(1.89E28*a**-1.55)/(4*pi*(a_planet*au2m*m2cm)**2.) < 6.3E-4*f for a,f in zip(age,flux)]
      tau_sat = age[argmax(temp)]
if diags:
      print("Sanz-Forcada et al. (2011) saturation time is "+str(round(tau_sat,4))+" Gyr.")
      if do_emp_sat and not do_emp_scale:
            print("        ---> Using empirically-fit 'saturation' time to match saturated/subsaturated fluxes from eq. 5")
      elif do_emp_scale:
            print("        ---> Using instead empirical fit for M dwarfs from Peacock et al. (2020)")

#X-ray (5-100 A) + EUV (100-920 A) = XUV (5-920 A)
xray_flux = []; euv_flux = []; xuv_flux = []
xray_scale = []; euv_scale = []; saturation = []
#Peacock et al. (2020) reports F_x/F_J and F_euv/F_J, which we scale with F_J_bol:
#     F_J/F_bol = 2.5119^(m_bol - m_J), where m_bol = 5log(d_star/10 pc) - 2.5log(L/L_sun) 
#                                         and m_J is stated in TOI-1266 discovery paper
if do_emp_scale:
      F_J_bol = 2.5119**(5*log10(d_star/10)-2.5*log10(l_star) - J_mag_star)
      if diags:
            print(f"F_J/F_bol = {F_J_bol:8.2e}")
for i in range(len(age)):
      if do_emp_scale: #empirical scaling from Peacock et al. (2020) - HAZMAT VI
            euv_lum = luminosity[i]*l_Sun*min(8.6E-3,2.38*(age[i]*1E3)**-1)*F_J_bol
            xray_lum = luminosity[i]*l_Sun*min(5.3E-3,10.78*(age[i]*1e3)**-1.36)*F_J_bol
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

#Handy lambda functions for dependent variables
f_roche = lambda rp,mp,ms,a : a*au2m*((mp*m_Earth)/(3*ms*m_Sun))**(1/3)/(rp*r_Earth)      #in planet radii
f_ktide = lambda ro,rp : (1 - (3/(2*(ro/rp))) + (1/(2*(ro/rp)**3)))                       #dimensionless
f_T_eq  = lambda f,a : ((f*(1-a))/(4*sigma))**0.25
f_grav  = lambda mp,rp : (G*mp*m_Earth)/(rp*r_Earth)**2.
f_mu    = lambda x,m : sum(x[i]*m[i] for i in range(len(x)))
f_H_below=lambda t,mu,g : kb*t/(mu*m_H*g)
f_rho_P = lambda pressure,temperature,mu : [mu*m_H*p/(kb*t) for p,t in zip(pressure,temperature)]
f_xuvcr = lambda mp,rp,eff: 180*mp**2*rp**-3*(eff/0.3)**-1 #Critical XUV flux for O escape from Luger & Barnes (2018)

#define initial variable states
c_f_t = [core_frac]; c_d_t = [core_den]
e_f_t = [envelope_frac]; e_comp_t = [envelope_comp]; num_comps = len(envelope_comp)
r_p_t = [r_planet]
m_p_t = [m_planet]
grav_t = [f_grav(m_planet,r_planet)] #Gravity [m/s2] at the 20 mbar level assumed to be for m_planet @ r_planet
roche_t = [f_roche(r_p_t[-1],m_p_t[-1],m_star,a_planet)]
ktide_t = [f_ktide(roche_t[-1],r_p_t[-1])]
mu_t = [f_mu(e_comp_t[-1],envelope_compm)]
vmr_t = [[e_comp_t[-1][i]*envelope_compm[i]/mu_t[-1] for i in range(len(envelope_compm))]]
epsilon_t = [f_mu(e_comp_t[-1],efficiencies)] #Re-use f_mu to get mean epsilon for composition
escape_flux_t = []; mescape_flux_t = []
crossover_mass_t = []; mass_change = 0.

if mode == 'forward':
      start = [0]
      end = [len(age)]
      step = [1]
      print("Forward mode currently not implemented. Try again tomorrow.")
elif mode == 'inverse':
      start = [current_age_ind,current_age_ind+1]
      end = [0,len(age)]
      step = [-1,1]
      print("Solving for initial planet mass from starting mass, radius, and composition in inverse mode.")
else:
      exit("Unrecognized mass/radius/composition mode.")

#Have to define this in here - T_eq depends on when we start.
T_eq_t = [f_T_eq(flux[start[0]]*ergcm2s2Wm2,albedo)]
      
if diags:
      print(f"Starting from {round(age[start[0]],3)} Gyr, L/L_sol = {luminosity[start[0]]:6.4f}, XUV/L_bol = {saturation[start[0]]:8.2e}")
      print(f"Critical XUV flux (O loss) = {f_xuvcr(m_planet,r_planet,epsilon_t[-1]):8.2e} erg/cm2/s; F_XUV = {xuv_flux[start[0]]:8.2e} erg/cm2/s")
      print('Starting values for some variables:')
      print('Roche = ',roche_t,' Ktide = ',ktide_t,' mu = ',mu_t,' T_eq = ',T_eq_t,' epsilon = ',epsilon_t)
      print('core_frac = ',c_f_t,' envelope_frac = ',e_f_t)
      print(f' envelope_comp: {envelope_species} = {envelope_comp}')
      
if estimate_atmosphere:
      atm_est(grav_t[-1],T_eq_t[-1],e_f_t[-1]*m_p_t[-1]*m_Earth,vmr_t[-1],mu_t[-1],diags,plotdir)

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
                  dts = age[i]*1E9*yr2s
            H_below = f_H_below(T_eq_t[-1],mu_t[-1],grav_t[-1])         #meters
            p_xuvbase = p_xuvcofac*grav_t[-1]                           #Pascal
            rbase = r_p_t[-1]*r_Earth + H_below*log(p_photo/p_xuvbase)  #meters
            #This would be a good place to call a function for determining which escape regime you're in, a la Owen & Alvarez (2016)
            if euv_flux[i] > 1E4:
                  #calculate escape rate with radiation-recombination limit from Lopez et al. (2017; MNRAS)
                  mflux = 7.11E4*(euv_flux[i])**0.5*r_p_t[-1]**1.5
            else:
                  #calculate energy-limited escape rate following Lopez et al. (2017; MNRAS)
                  mflux = epsilon_t[-1]*pi*xuv_flux[i]*ergcm2s2Wm2*rbase**3./(G*m_p_t[-1]*m_Earth*ktide_t[-1])
            eflux = mflux/(m_H*4*pi*(r_p_t[-1]*r_Earth)**2.) #assume the escape flux is H to start, then iterate through comp_mass list
            ref_H_flux = epsilon_t[-1]*xuv_flux[i]*ergcm2s2Wm2*r_p_t[-1]*r_Earth/(4*G*m_p_t[-1]*m_Earth*ktide_t[-1]*m_H)
            if finediags:
                  print(f"L&B15 F_H_ref = {ref_H_flux:8.2e} vs. escape flux = {eflux:8.2e} [molec/m2/s]")
            crossover_mass_t.append(crossover_mass(T_eq_t[-1],eflux,grav_t[-1],2*e_comp_t[-1][0],mass=envelope_compm[0]/2,mu=mu_t[-1]))
            for m in range(num_comps):
                  if crossover_mass_t[-1] > envelope_compm[m]:
                        break
            if finediags:
                  print(f"Component {m} has mass {envelope_compm[m]:4.3f} which is greater than "+
                        f"crossover mass = {crossover_mass_t[-1]:4.3f} for flux = {eflux:8.2e} molec/m2/s")
            mescape_flux_t.append(mflux)
            escape_flux_t.append(eflux)
            if finediags:
                  print(f'Escape flux @ {age[i]:5.3f} Gyr = {escape_flux_t[-1]:8.2e} molec/m2/s ({mflux:8.2e} kg/s)') 
            
            #update variables
            #Append new values as we move backwards in time, then 'flip' variables
            mass_H = -1*step[j]*mflux*dts/m_Earth #Earth masses
            new_mass = [(j==0)*mass_H + e_comp_t[-1][j]*e_f_t[-1]*m_p_t[-1] for j in range(num_comps)]
            #if any([m<0 for m in new_mass]):
            #      exit("Mass of one of the reservoirs is negative."+str(new_mass))
            new_e_frac = [m/sum(new_mass) for m in new_mass]
            mass_change = mass_change+mflux*dts
            m_p_t.append(m_p_t[-1]+mass_H)
            c_f_t.append(core_frac*m_planet/m_p_t[-1])
            e_f_t.append(1.-c_f_t[-1])
            e_comp_t.append(new_e_frac)
            r_p_t.append(r_p_t[-1])
            grav_t.append(f_grav(m_p_t[-1],r_p_t[-1]))
            roche_t.append(f_roche(r_p_t[-1],m_p_t[-1],m_star,a_planet))
            ktide_t.append(f_ktide(roche_t[-1],r_p_t[-1]))
            mu_t.append(f_mu(e_comp_t[-1],envelope_compm))
            vmr_t.append([e_comp_t[-1][i]*envelope_compm[i]/mu_t[-1] for i in range(len(envelope_compm))])
            epsilon_t.append(f_mu(e_comp_t[-1],efficiencies)) #Re-use mubar calculator to get mean epsilon
            T_eq_t.append(f_T_eq(flux[i]*ergcm2s2Wm2,albedo))
            if finediags:
                  print(f'mass added = {mass_H} Earth masses, dts = {dts}, efficiency = {epsilon_t[-1]}, envelope composition = {e_comp_t[-1]}')
            if i%10 == 5 and finediags:
                  exit()
      
      if mode == 'inverse' and j == 0:
            #Flip variables so that they're aligned with *age* variable
            c_f_t.reverse(); c_d_t.reverse()
            e_f_t.reverse(); e_comp_t.reverse()
            r_p_t.reverse(); m_p_t.reverse()
            mu_t.reverse(); T_eq_t.reverse()
            roche_t.reverse(); ktide_t.reverse()
            epsilon_t.reverse(); crossover_mass_t.reverse()
            mescape_flux_t.reverse()
            escape_flux_t.reverse()


e_cmp_t = [[x*a*b for x in e_comp_t[e_f_t.index(a)]] for a,b in zip(e_f_t,m_p_t)]
mescape_flux_t.insert(0,mescape_flux_t[0])
escape_flux_t.insert(0,escape_flux_t[0])
crossover_mass_t.insert(0,crossover_mass_t[0])

if diags:
      print(f"{max(m_p_t):6.4f}-->{min(m_p_t):6.4f} m_Earth; total mass change = {mass_change/m_Earth:6.4f} m_Earth")
      print(f"Present-day escape rate = {mescape_flux_t[age.index(current_age)]:8.2e} kg/s")
      print(f"                        = {escape_flux_t[age.index(current_age)]:8.2e} molec/m2/s")
      print(f"                        = {escape_flux_t[age.index(current_age)]/1e4:8.2e} molec/cm2/s")
      #print(c_f_t[0],e_f_t[0],c_f_t[-1],e_f_t[-1])
      #print(e_comp_t[0],e_comp_t[-1])
      #print(e_cmp_t[0],e_cmp_t[-1])
      #print(min(crossover_mass_t),sum(crossover_mass_t)/len(crossover_mass_t),max(crossover_mass_t),' Crossover mass evolution')

if plots:
      fig1,ax1 = plt.subplots()
      if diags:
            ax1.plot(a1,l1,'r--') #raw Baraffe data
            if interp_lum:
                  ax1.plot(a2,l2,'r--') #other Baraffe data
                  ax1.plot([1.1*max(age),1.1*max(age)],[0.98*f_lum_CW18(m_star),1.02*f_lum_CW18(m_star)],color='gray')
      ax1.plot(age,luminosity,'k')
      #include saturation timescale
      ax1.plot(tau_sat,luminosity[age.index(tau_sat)],'bo')
      ax1.text(tau_sat,1.1*luminosity[age.index(tau_sat)],r'${\tau}_{sat}$',color='b',ha='left')
      #include the notional start of the main sequence lifetime for star
      ax1.plot(ms_start,luminosity[age.index(ms_start)],'ro',markerfacecolor='None')
      ax1.text(ms_start,1.1*luminosity[age.index(ms_start)],'MS',color='r',ha='right')
      #include the notional age of the system
      ax1.plot(current_age,luminosity[current_age_ind],'x',color='gray')
      ax1.text(current_age,1.1*luminosity[current_age_ind],'Today',color='gray',ha='center')
      ax1.set_xscale('log'); ax1.set_yscale('log')
      ax1.set_xlabel("Age [Gyr]"); ax1.set_ylabel("Luminosity [L/L$_{\odot}$]")
      ax1.set_title("Baraffe et al. (2015) luminosity evolution for "+str(match)+" M$_{\odot}$ star\n"
            +"Raw Baraffe data(r--) and interpolation(k)")
      plt.show()
      if save_plots:
            fig1.savefig(plotdir+'luminosity_MS_tau_i_vs_time.png',dpi=200)

if plots:
      fig2,ax2 = plt.subplots()
      ax2.plot(age,m_p_t,age,[i*j for i,j in zip(m_p_t,e_f_t)],age,[i*j for i,j in zip(m_p_t,c_f_t)])
      ax2.set_xscale('log');ax2.set_yscale('log')
      ax2.set_xlabel("Age [Gyr]"); ax2.set_ylabel(r"Planet Mass [$M_{\oplus}$]")
      ax2.plot(current_age,m_p_t[current_age_ind],'x',color='gray')
      ax2.text(current_age,1.02*m_p_t[current_age_ind],'Today',color='gray',ha='center')
      ax2.text(age[-1],m_p_t[-1]*e_f_t[-1],r"M$_{atm}$")
      ax2.text(age[-1],m_p_t[-1]*c_f_t[-1],r"M$_{core}$")
      #for i in range(num_comps):
      #      ax2.text(age[-1],e_cmp_t[-1][i],envelope_species[i])
      plt.show()
      if save_plots:
            fig2.savefig(plotdir+'Planet_mass_vs_time.png',dpi=200)

if plots:
      fig3,ax3 = plt.subplots()
      ax3.plot(age,saturation,'k--')
      ax3.plot(age,xray_scale,'r')
      ax3.plot(age,euv_scale,'b')
      ax3.set_xscale('log'); ax3.set_yscale('log')
      ax3.set_xlabel("Age [Gyr]"); ax3.set_ylabel(r"X-ray(r)/EUV(b)/XUV(k--) sat. [L$_{\lambda}$/L$_{bol}$]")
      ax3.set_title("Baraffe et al. (2015) X-ray, EUV, and XUV evolution "+str(match)+" M$_{\odot}$ star")
      plt.show()
      if save_plots:
            fig3.savefig(plotdir+'UV_scaling_vs_time.png',dpi=200)

if plots:
      fig5,ax5 = plt.subplots()
      ax5.plot(age,xuv_flux,'k--')
      ax5.plot(age,xray_flux,'r')
      ax5.plot(age,euv_flux,'b')
      ax5.plot([min(age),max(age)],[1e4,1e4],color='gray',linestyle='-.')
      print(f"XUV flux drops below 1E4 erg/cm2/s at {age[find_nearest(xuv_flux,1E4)]} Gyr.")
      ax5.set_xscale('log'); ax5.set_yscale('log')
      ax5.set_xlabel("Age [Gyr]"); ax5.set_ylabel(r"X-ray(r)/EUV(b)/XUV(k--) flux [ergs/cm2/s]")
      ax5.set_title("Baraffe et al. (2015) X-ray, EUV, and XUV evolution "+str(match)+" M$_{\odot}$ star")
      plt.show()
      if save_plots:
            fig5.savefig(plotdir+'UV_flux_vs_time.png',dpi=200)

if plots:
      fig6,ax6 = plt.subplots()
      ax6.plot(age,crossover_mass_t)#test_mass)#[e/1e4 for e in escape_flux_t])#test_mass)#mescape_flux_t)
      ax6.set_xscale("log"); ax6.set_yscale('log')
      ax6.set_xlabel("Age [Gyr]"); ax6.set_ylabel('Crossover mass [amu]')#'Escape flux [molec/cm2/s]')#'Change in mass [kg]')#r"Escape rate [kg/s]")
      crit_mass = find_nearest(crossover_mass_t,18.)
      ax6.plot(age[crit_mass-1:crit_mass+1],crossover_mass_t[crit_mass-1:crit_mass+1],'k-')
      ax6.text(age[crit_mass],1.1*crossover_mass_t[crit_mass],"Water")
      plt.show()
      if save_plots:
            fig6.savefig(plotdir+'crossover_mass_vs_time.png',dpi=200)#'escape_rate_vs_time.png',dpi=200)#'mass_loss_vs_time.png',dpi=200)
