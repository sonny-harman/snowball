from math import pi,exp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#case-specific variables, constants, and functions
from planet import p_photo,envelope_species,envelope_compm
from planet import r_planet,m_planet
from constants import G,kb,m_H,Rconst,sigma
from constants import r_Earth,m_Earth
from modules import read_thermo,read_pt_profile

a_H2 = [4.966884120E+08,-3.147547149E+05,7.984121880E+01,-8.414789210E-03,4.753248350E-07,-1.371873492E-11,1.605461756E-16,2.488433516E+06,-6.695728110E+02]

#local variable definitions
H_frac = 0.01
match_rc = True

def atm_est(g,t_eq,envelope_mass,r_core,vmr,mu,diags,plotdir):
      #Read in thermochemical data for calculating c_p of species
      coeffs = {}
      for species in envelope_species:
            if diags:
                  print('Retrieving thermochemical data for '+species)
            coeffs[species] = read_thermo(species,diags)
      #c_pm = [J/K/kg] = c_p / mu = [J/K/mol] / [kg/mol]
      if mu == 2.:
            f_c_pm = lambda n,T : ( ((T<1000)*(sum([coeffs[n][0][i+2]*T**i for i in range(-2,5)])) 
                                  + (T>=1000 and T<6000)*(sum([coeffs[n][1][i+2]*T**i for i in range(-2,5)])) 
                                  + (T>=6000 and T<2E4)*(sum([a_H2[i+2]*T**i for i in range(-2,5)])))
                                  / (envelope_compm[envelope_species.index(n)]/1000.) )
            T_max = 2E4
      else:
            f_c_pm = lambda n,T : ( ((T<1000)*(sum([coeffs[n][0][i+2]*T**i for i in range(-2,5)])) 
                                  + (T>=1000 and T<6000)*(sum([coeffs[n][1][i+2]*T**i for i in range(-2,5)])))
                                  / (envelope_compm[envelope_species.index(n)]/1000.) )
            T_max = 6E3
      f_c_pmtot= lambda x,T : sum([f_c_pm(envelope_species[i],T)*x[i] for i in range(len(x))])
      f_grav   = lambda mp,rp : (G*mp*m_Earth)/(rp*r_Earth)**2.
      f_H_below= lambda t,mu,g : kb*t/(mu*m_H*g)

      #calculate the extent of the atmosphere, assuming a dry adiabat
      K_zed = []
      p_temp = [p_photo]; alt_temp = [0.]; T_temp = [t_eq]
      g_temp = g
      m_above = p_photo*4*pi*(r_planet*r_Earth)**2/g_temp
      rk_p,rk_alt,rk_t = read_pt_profile()
      if match_rc:
            f_temp_interp = interp1d(rk_p,rk_t,kind='cubic')
            T_temp[-1] = f_temp_interp(p_photo)
            if diags:
                  print(f"Matching data in P/T profile at {p_photo:8.2e} bars, where T = {T_temp[-1]:6.1f}K")
      if diags:
            print(f"Mass above {p_photo/1E5:6.0e} bars = {m_above:10.2e} kg; total envelope is {envelope_mass:10.2e}; T_eq = {T_temp[-1]:7.2f}K")
      while m_above < envelope_mass and T_max > T_temp[-1]:
            rp_temp = r_planet*r_Earth-alt_temp[-1]
            mp_temp = m_planet*m_Earth-m_above
            g_temp = f_grav(mp_temp/m_Earth,rp_temp/r_Earth)
            m_above = p_temp[-1]*4*pi*rp_temp**2/g_temp
            H_temp = f_H_below(T_temp[-1],mu,g_temp)
            c_temp = f_c_pmtot(vmr,T_temp[-1])
            Gamma_temp = g_temp/c_temp
      #now we can update to a new pressure/altitude level
            alt_temp.append(alt_temp[-1] + H_temp*H_frac)
            T_temp.append(T_temp[-1] + Gamma_temp*H_temp*H_frac)
            #T_temp.append(min(5999,T_temp[-1] + Gamma_temp*H_temp*H_frac)) #thermochemical data is only good to 6,000 K
            p_temp.append(p_temp[-1]*exp(H_frac))
            #Calculating the Kzz from Gao et al., 2018 which lays out Ackerman and Marley (2001)'s methodology - produces Kzz ~100 m2/s (1E6 cm2/s), which is low
            #L_mix = max(0.1*H_temp,H_temp*(T_temp[-1]-T_temp[-2])/(f_temp_interp(p_temp[-1])-f_temp_interp(p_temp[-2])))
            #Kzz = (H_temp/3)*(L_mix/H_temp)**(4/3)*(Rconst**2*sigma*T_temp[-1]**4/(mu**2*p_temp[-1]*c_temp*T_temp[-1]))**(1/3)
            #K_zed.append(Kzz)
            if r_planet*r_Earth - alt_temp[-1] < 0.:
                  print('Run out of room - no more planet!')
                  print(rp_temp,alt_temp[-1],T_temp[-1])
                  break
      if diags:
            print(f"Mass above {p_temp[-1]/1E5:8.2f} bars = {m_above:10.2e} kg, or {m_above/m_Earth:7.4f} Earth masses")
            print(f"g({p_temp[-1]/1E5:8.2f} bar) = {g_temp:6.2f} m/s2, H_atm = {H_temp/1.E3:7.2f} km, delta_r/r_p = {100*alt_temp[-1]/(r_Earth*r_planet):5.2f}% of r_p, T = {T_temp[-1]:7.2f}")
            print(f"Dry adiabatic lapse rate @ {p_temp[-1]/1E5:8.2f} bars = {Gamma_temp*1000:5.2f} K/km")
            fig4,ax4 = plt.subplots()
            #ax4.semilogy(T_temp,[p/1E5 for p in p_temp],rk_t,[p/1E5 for p in rk_p])
            #ax4.set_xlabel('Pressure [bar]'); ax4.set_ylabel('Temperature [K]')
            #ax4.invert_yaxis()
            ax4.loglog(T_temp,p_temp,rk_t,rk_p)
            ax4.set_ylabel('Pressure [Pa]'); ax4.set_xlabel('Temperature [K]')
            ax4.set_xlim(200,10000); ax4.set_ylim(100,1e12)
            ##rho1 = f_rho_P(p_temp,T_temp,mu_t[-1])
            ##rho2 = f_rho_P(rk_p,rk_t,mu_t[-1])
            ###ax4.plot(rho1,[p/1e5 for p in p_temp],rho2,[p/1e5 for p in rk_p])
            ###ax4.set_xlabel('Density [kg/m3]'); ax4.set_ylabel('Pressure [bar]')
            ##ax4.plot(rho1,T_temp,rho2,rk_t)
            ##ax4.set_xlabel('Density [kg/m3]'); ax4.set_ylabel('Temperature [K]')
            ####ax4.set_xlim(0, 800); ax4.set_ylim(0,1000)
            #ax4.loglog([p/1E8 for p in p_temp],T_temp,[p/1e8 for p in rk_p],rk_t)
            #ax4.set_xlabel('Pressure [kbar]'); ax4.set_ylabel('Temperature [K]')
            ax4.set_title(f"{m_above/m_Earth:6.2f} Earth-mass atmosphere\n {[s+'='+str(round(v,4)) for s,v in zip(envelope_species,vmr)]} [VMR]")
            fig4.savefig(plotdir+'atmospheric_estimate.png',dpi=200)
            plt.show()
      
            #fig5,ax5 = plt.subplots()
            #ax5.loglog(K_zed,p_temp[1:])
            #ax5.set_ylabel('Pressure [Pa]'); ax5.set_xlabel('Kzz [/m2/s]')
            #plt.show()
      return


