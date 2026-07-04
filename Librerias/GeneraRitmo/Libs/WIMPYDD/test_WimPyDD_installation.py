import WimPyDD as WD
import numpy as np
import matplotlib.pyplot as pl

### if WimPyDD is downloaded using the clone command the modification time of all files is lost
### and set to the cloning time.
### As a consequence  when the routine load_response_functions is called  WimPyDD issues a warning because
### the files with the integrated response functions do not appear to have been created after 
### the input files containing the energy resolution, quenching, binning, exposure etc.
### calling the routine update_response_functions_time_stamps() solves this issue.
WD.update_response_functions_time_stamps()

pl.clf()
xe=WD.Xe

### defines an array with recoil energy values
er=np.linspace(0,1400,1000) # in keV

## loop over isotopes of Xe (the attributes xe.mass,xe.isotopes,xe.func_w are all arrays with one entry per isotope)
for m,isotope,func_w in zip (xe.mass,xe.isotopes,xe.func_w):
    q_vec=np.sqrt(2*m*er*1e-6)# gets corresponding transferred monentum in GeV (required input for func_w functions) 
    nn=WD.nuclear_current['Phi_prime_prime']# check WD.nuclear_currents.keys()
    tau=0
    tau_prime=0
    
    func_w_phi_prime_prime=[func_w(q)[nn,tau,tau_prime] for q in q_vec]
    pl.plot(er,func_w_phi_prime_prime,label=isotope)

pl.title('$W^{00}_{\Phi^{\prime\prime}}$')
pl.yscale('log')
pl.legend()
pl.show()

pl.clf()
nai=WD.DAMA_LIBRA_2019.target
na,i=nai.element
#N.B. using the DAMA_LIBRA target sodium and iodine have the quenching attribute:
# 'quenching' in dir(na) -> True
# 'quenching' in dir(i) -> True
# the routine diff_rate interprets the argument er as electron-equivalent energy
# E_ee=Q(E_R)*E_R with Q the quenching and E_R the nuclear recoil energy.

# On the other hand:
# 'quenching' in dir(WD.Na) -> False
# 'quenching' in dir(WD.I) -> False
# i.e. the pre-defined targets do not have quenching information.
# if used in diff_rate the argument er is interpreted as nuclear recoil energy. 

#Define Hamiltonian in terms of two parameters: the effective scale M and r=cn/cp
# (ratio  of WIMP-neutron over WIMP-proton couplings) 
wc={1: lambda M, r=1 : [(1.+r)/M**2, (1-r)/M**2] }
SI_M=WD.eft_hamiltonian('SI_M',wc)

er_vec=np.linspace(0.1,20,100) 
mchi=20. 
vmin,delta_eta0=WD.streamed_halo_function()# Maxwellian with standard parameters
diff_rate=[WD.diff_rate(nai, SI_M, mchi, er, vmin, delta_eta0,M=1e3) for er in er_vec]
# N.B the r parameter is not passed, default value r=1 used. 
# N.B.2 also the argument exposure is not passed, so the rate is normalized to events/kg/day/keV. 
pl.plot(er_vec,diff_rate)

# loads response functions (j_chi=0.5 spin default value)
# from the directory WimPyDD/Experiments/DAMA_LIBRA_2019/Response_functions/spin_1_2/
# since SI_M.coeff_squared_list -> [(1, 1)]
# the file c_1_c_1.npy is loaded.

WD.load_response_functions(WD.DAMA_LIBRA_2019,SI_M,verbose=False)
# before calling load_response_functions the dictionary attribute
# WD.DAMA_LIBRA_2019.response_functions is empty. After the call:
# [(hamiltonian.name,spin) for hamiltonian,spin in WD.DAMA_LIBRA_2019.response_functions.keys()] -> [('SI_M', 0.5)], i.e. the response functions for the SI_M hamiltonian and spin 1/2 are loaded.

e_prime,s0=WD.wimp_dd_rate(WD.DAMA_LIBRA_2019, SI_M, vmin, delta_eta0, mchi,M=1e3)
#N.B the output is in events/kg/day/keV because in the WimPyDD/Experiments/DAMAthe exposure is set
# to 1/delta_e, with delta_e the with of each bin:
#WD.DAMA_LIBRA_2019.exposure  -> [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.125]

pl.step(e_prime,s0)

pl.xlabel('E$^{\prime}$ or $E_{ee}$ (keV)')
pl.ylabel(' $dR/dE_R$ (events/kg/day/keV)') 

pl.yscale('log')

pl.show()


pl.clf()
vmin,delta_eta1=WD.streamed_halo_function(yearly_modulation=True)
# delta_eta1 for yearly modulation calculated setting yearly_modulation=True
#

# define Wilson coefficient c9 in terms of the reference cross section
# sigma_ref=c9^2*mu_chi_n^2/pi, with mu_chi_n the WIMP-nucleon reduced mass, # and of r=c9^n/c9^p (WIMP-neutron over WIMP-proton ratio)
# sigma_ref in cm^2 and c9 in GeV^-2.

def c9_tau(mchi,cross_section,r):
    hbarc2=0.389e-27
    mn=0.931
    mu=mchi*mn/(mchi+mn)
    cp=(np.pi*cross_section/(mu**2*hbarc2))**(1./2.)
    return cp*np.array([1.+r,1.-r])

# instantiates the effective Hamiltonian containing only O9.
c9=WD.eft_hamiltonian('c9',{9: c9_tau})

# best-fit values from Table 1 of rXiv:1804.07528.
cross_section=8.29e-33
mchi=9.3
r=4.36

# load the response functions (j_chi=0.5 is the spin default value) 
WD.load_response_functions(WD.DAMA_LIBRA_2019,c9,verbose=False)

e,sm=WD.wimp_dd_rate(WD.DAMA_LIBRA_2019, c9, vmin, delta_eta1,cross_section=cross_section, mchi=mchi , r=r)

pl.plot(e,sm,'-k',linewidth=3)

# modulation amplitudes measured by DAMA-LIBRA (from Fig. 11 of NUCL. PHYS. AT. ENERGY 19 (2018) 307-325)
experimental_s1=np.array([0.0232    , 0.0164    , 0.0178    , 0.019     , 0.0178,0.0109, 0.011, 0.004, 0.0065, 0.0066,0.0009,0.0016,0.0007,0.0016, 0.00039922])

errors_on_s1=np.array([0.0052, 0.0043, 0.0028    , 0.0029    , 0.0028,0.0025    , 0.0022    , 0.002     , 0.002, 0.0019,0.0018    , 0.0018    , 0.0018    , 0.0018    , 0.00046397])

errors_on_energy=np.array([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,0.25, 0.25, 0.25, 4.  ])
pl.errorbar(e,experimental_s1,yerr=errors_on_s1,xerr=errors_on_energy,fmt='none',ecolor='red')

pl.xlabel('E$^{\prime}$ (keV)')
pl.ylabel('$S^{1}_{[E_1^{\prime},E_2^{\prime}]}$ (events/kg/day/keV)')

pl.show()

pl.clf()
vmin,delta_eta1=WD.streamed_halo_function(yearly_modulation=True)
# delta_eta1 for yearly modulation calculated setting yearly_modulation=True
#
# defines c_tau (in GeV^-2) in terms of the reference cross section
# sigma_ref=|c_vec|^2*mu_chi_n^2/pi (in cm^2)
# where c_vec=[c4^1,c4^2,c5^1,c5^2,c6^1,c6^2]
# and mu_chi_n is the WIMP-nucleon reduced mass.
def c_tau(mchi,sigma_ref,c_tilde):
    hbarc2=0.389e-27
    mn=0.931
    mu=mchi*mn/(mchi+mn)
    c_abs=(mu**2*hbarc2/(np.pi*sigma_ref))**(-0.5)
    return c_abs*c_tilde

# c_tilde=c_vec/|c_vec| is the unit vector pointing in the direction of c_vec
# that minimizes the tension with the DAMA modulation amplitudes
c_tilde=np.array([-0.00139183, -0.00145925,-0.00318353, -0.01660634, 0.69199713, 0.72169939])

# values of WIMP mass, mass splitting and sigma_ref that minimize the tension with DAMA
# modulation amplitudes
mchi=11.64102564102564
delta=23.73913043478261
sigma_ref=4.684124873205523e-28
# calculates c_vec in terms of sigma_ref, mchi,c_tilde
c_vec=c_tau(mchi,sigma_ref,c_tilde)

# defines a Hamiltonian that maps mchi, sigma_ref and c_tilde into the couplings c_4^tau, c_5^tau, c_6^tau

wc={4: lambda mchi,s,c: c_tau(mchi,s,c)[0:2],5: lambda mchi,s,c: c_tau(mchi,s,c)[2:4],6: lambda mchi,s,c: c_tau(mchi,s,c)[4:6]}
c_456=WD.eft_hamiltonian('c_456',wc)

# loads the response functions (j_chi=0.5 is default valu for spin)
WD.load_response_functions(WD.DAMA_LIBRA_2019,c_456,verbose=False)

# calculates the modulation amplitudes in terms of mchi, sigma_ref and c_tilde
e_vec,dama_th=WD.wimp_dd_rate(WD.DAMA_LIBRA_2019,c_456,vmin,delta_eta1,mchi,delta=delta,s=sigma_ref,c=c_tilde)
pl.plot(e_vec,dama_th,'-k',linewidth=1)


# modulation amplitudes measured by DAMA-LIBRA (from Fig. 11 of NUCL. PHYS. AT. ENERGY 19 (2018) 307-325)
experimental_s1=np.array([0.0232    , 0.0164    , 0.0178    , 0.019     , 0.0178,0.0109, 0.011, 0.004, 0.0065, 0.0066,0.0009,0.0016,0.0007,0.0016, 0.00039922])

errors_on_s1=np.array([0.0052, 0.0043, 0.0028    , 0.0029    , 0.0028,0.0025    , 0.0022    , 0.002     , 0.002, 0.0019,0.0018    , 0.0018    , 0.0018    , 0.0018    , 0.00046397])

errors_on_energy=np.array([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,0.25, 0.25, 0.25, 4.  ])
pl.errorbar(e,experimental_s1,yerr=errors_on_s1,xerr=errors_on_energy,fmt='none',ecolor='red')

pl.xlabel('E$^{\prime}$ (keV)')
pl.ylabel('$S^{1}_{[E_1^{\prime},E_2^{\prime}]}$ (events/kg/day/keV)')

pl.show()


pl.clf()
vmin, delta_eta0=WD.streamed_halo_function()


# defines c_1^tau=[c_1^0, c_1^1] in terms of the WIMP-nucleon cross section sigma_N=(c_N)^2*mu_chi_N^2/pi, with c_N=c_1^p=c_1^n, c_1^0=c_1^p+c_1^n,
# c_1^1=c_1^p-c_1^n (for r=c_1^n/c_1^p=1)
def c_tau_SI(sigma_N,mchi,cn_over_cp=1):
    hbarc2=0.389e-27 #(hbar*c)^2 in GeV^2 * cm^2
    mn=0.931
    mu=mchi*mn/(mchi+mn)
    return np.sqrt(np.pi*sigma_N/hbarc2)/mu*np.array([1+cn_over_cp,1-cn_over_cp])

SI=WD.eft_hamiltonian('Spin-independent',{1: c_tau_SI})

WD.load_response_functions(WD.XENON_1T_2018,SI,verbose=False)

# uses built-in XENON_1T_2018 experiment 
mchi,sigma=WD.mchi_vs_exclusion(WD.XENON_1T_2018, SI,vmin, delta_eta0, n_points=20,mchi_scale='log')
pl.plot(mchi,sigma,':r',label='WimPyDD')


# published exclusion plot, from arXiv:1805.12562
mchi_xenon1t=np.array([  6.38084285,   6.82568656,   7.38401692,   8.16949332,
         9.14061906,  10.1129544 ,  11.44291272,  13.39150131,
        16.76448605,  19.61927374,  23.2195425 ,  27.48048481,
        31.80087168,  36.3894599 ,  45.04623899,  54.52370581,
        74.67422632, 102.2718466 , 160.27952577, 262.73153013,
       376.36545636, 539.14715404, 705.96063944])

sigma_xenon1t=np.array([1.34339933e-44, 7.44380301e-45, 3.66524124e-45, 1.80472177e-45,
       9.15247311e-46, 5.37983840e-46, 2.89426612e-46, 1.65176674e-46,
       8.37677640e-47, 6.23550734e-47, 4.78065253e-47, 4.24820170e-47,
       4.24820170e-47, 4.37547938e-47, 4.92388263e-47, 5.54102033e-47,
       7.22727132e-47, 9.42668455e-47, 1.42510267e-46, 2.28546386e-46,
       3.25702066e-46, 4.64158883e-46, 5.87801607e-46])

pl.plot(mchi_xenon1t,sigma_xenon1t,'-k',label='Published')
pl.title('SI, XENON1T exclusion plot')
pl.xscale('log')
pl.yscale('log')
pl.xlabel('$m_\chi$ (GeV)')
pl.ylabel('$\sigma_p$ (cm$^2$)')
pl.legend()
pl.show()



pl.clf()
vmin, delta_eta0=WD.streamed_halo_function()

pl.clf()

# defines the effective Hamiltonian for spin-independent interaction parameterizing the Wilson coefficients
# interms of the effective scale M: c_1^tau=[c_1,c_1^1]=[c_1^p+c_1^n,c_1^p+c_1^n] with c_1^p=c_1^n=1/M^2
SI_M=WD.eft_hamiltonian('SI',{1: lambda M:1./M**2*np.array([2,0])})

# loads response functions for built-in experiment WD.XENON_1T_2018 (j_chi=0.5 defalult spin value)
WD.load_response_functions(WD.XENON_1T_2018,SI_M,verbose=False)

# calling mchi_vs_exclusion without passing M corresponds to setting M=1 (holds for any parameter without a default value
# the rate is proportional to 1/M^4, so the exclusion plot is on the same quantity.
mchi,one_over_M4=WD.mchi_vs_exclusion(WD.XENON_1T_2018, SI_M, vmin,delta_eta0,n_points=20,mchi_scale='log')

# plots M=(1/M^4)^(-1/4)
pl.plot(mchi,one_over_M4**(-1./4.))
pl.xlabel('$m_\chi$ (GeV)')
pl.ylabel('$M$ (GeV)')
pl.xscale('log')
pl.yscale('log')
pl.show()

##### defines a Hamiltonian containing the full set of operators in the base of Anand et al. and test
##### WD.load_response_functions on the three pre-defined experiment objects WD.XENON_1T_2018,
###### WD.DAMA_LIBRA_2019, WD.PICO60_2019

## the full set of operators 1,3,...,15
coupling_list=np.delete(np.arange(1,16),1)

# defines hamiltonian with "dummy" Wilson coefficients -   just to test correct loading of tables
wc={}
for coupling in coupling_list:
    wc.update({coupling:lambda :[1,1]})
hamiltonian=WD.eft_hamiltonian('all_couplings',wc)

WD.load_response_functions(WD.XENON_1T_2018,hamiltonian)
WD.load_response_functions(WD.DAMA_LIBRA_2019,hamiltonian)
WD.load_response_functions(WD.PICO60_2019,hamiltonian)

# plots the full set of nonvanishing response functions for XENON1T
WD.plot_response_functions(WD.XENON_1T_2018,hamiltonian)
pl.show()


# plots the full set of nonvanishing response functions for DAMA-LIBRA
WD.plot_response_functions(WD.DAMA_LIBRA_2019,hamiltonian)
pl.show()

# plots the full set of nonvanishing response functions for PICO60
WD.plot_response_functions(WD.PICO60_2019,hamiltonian)
pl.show()

print(hamiltonian)

# plots the H(v) response function for halo-independent calculations
vmin=np.linspace(1e-3,782,100)
pl.clf()
energy,h=WD.wimp_dd_rate(WD.XENON_1T_2018, WD.SI, vmin, 1/vmin, 100,sum_over_streams=False,sigma_p=9.066504208570347e-47)
pl.plot(vmin,h[0])
pl.text(500,15,r'$m_{\chi}$ = 100 GeV',fontsize=12)
pl.text(500,13,r'${\sigma}^p = 9.06 \times 10^{-47} cm^2$',fontsize=12)
pl.xlabel(r'$v_{min}$ (km/s)',size=12)
pl.ylabel('H(v) (events)',size=12)
pl.show()
