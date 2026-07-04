import WimPyDD as WD
import math
import numpy as np
import matplotlib.pyplot as pl
from timeit import default_timer as timer


### density
pl.clf()
pl.plot(WD.Sun.r_vec,WD.Sun.rho_tot,label='Total density')
pl.plot(WD.Sun.r_vec,WD.Sun.rho_i['1H'],label='1H')
pl.xlabel('r',fontsize=12)
pl.ylabel(r'$\rho$',fontsize=12)
pl.legend()
pl.show()


### response functions
WD.load_response_functions_capture(WD.Sun,WD.SI,verbose=False)

pl.clf()
WD.plot_response_functions_capture(WD.Sun,WD.SI)
pl.show()


### Sun following 1501.03729
ti=timer()

mv=246.2 ## GeV

# instantiates the effective Hamiltonian containing only O1.
model_SI=WD.eft_hamiltonian('SI',{1: lambda c0, c1: np.array([c0,c1])})
WD.load_response_functions_capture(WD.Sun,model_SI)

# Halo function
vmin, delta_eta=WD.streamed_halo_function()

mchi_vec=np.logspace(1, 3, 20) # in GeV

capture_vec=[]

print('Sun, spin-independent interaction:')

for n,mchi in enumerate(mchi_vec):
    s1=timer()
    capture=WD.wimp_capture(WD.Sun, model_SI, vmin=vmin, delta_eta=delta_eta, mchi=mchi, rho_loc=0.4, c0=1e-3/mv**2, c1=0,targets_list=[WD.H1])
    s2=timer()
    print('mchi=',mchi,'('+str((n+1)/20*100)+'%) completed in ',s2-s1,' s')
    capture_vec=np.append(capture_vec,capture)

pl.clf()
pl.plot(mchi_vec,capture_vec,color='darkblue',label='WimPyC')
pl.legend(frameon=False, fontsize=15)
pl.xscale('log')
pl.yscale('log')
pl.xlim(mchi_vec[0],mchi_vec[-1])
pl.tick_params(axis='both', which='major', labelsize=16)  # Increase major tick labels
pl.tick_params(axis='both', which='minor', labelsize=12)
pl.tight_layout(pad=2.2)
pl.xlabel(r'$m_\chi$ [GeV]', fontsize=15)
pl.ylabel(r'C [1/s]', fontsize=15)
pl.title(r'$C_' + str(1) + r'$, $^{1}$H', fontsize=15)
tf=timer()
print('Total time taken:',tf-ti)
pl.show()


### Sun long-range following 1305.0912
ti=timer()

def c_tau_SI_qm2(sigma_p,mchi,q,cn_over_cp=1):
    hbarc2=0.389e-27 #(hbar*c)^2 in GeV^2 * cm^2
    mn=0.931
    mu=mchi*mn/(mchi+mn)
    return np.sqrt(np.pi*sigma_p/hbarc2)/mu*np.array([(1+cn_over_cp)/q,(1-cn_over_cp)/q])
    

def c_tau_SI_qm4(sigma_p,mchi,q,cn_over_cp=1):
    hbarc2=0.389e-27 #(hbar*c)^2 in GeV^2 * cm^2
    mn=0.931
    mu=mchi*mn/(mchi+mn)
    return np.sqrt(np.pi*sigma_p/hbarc2)/mu*np.array([(1+cn_over_cp)/q**2,(1-cn_over_cp)/q**2])
    
# instantiates the effective Hamiltonian containing only O1/SI long-range.    
model_qm2=WD.eft_hamiltonian('Spin-independent(qm2)',{(1,'qm2'): c_tau_SI_qm2})
WD.load_response_functions_capture(WD.Sun,model_qm2)

mchidd_vec1=np.array([4.14979e+0,4.31711e+0,4.53578e+0,4.76553e+0,5.05665e+0,5.52700e+0,6.10109e+0,6.73482e+0,7.28888e+0,7.81095e+0,8.45354e+0,9.23986e+0,1.05065e+1,1.18294e+1,1.38558e+1,1.75645e+1,2.36260e+1,2.93637e+1,3.87242e+1,5.74985e+1,8.04599e+1,1.11483e+2,1.51446e+2,2.07777e+2,3.02473e+2,4.06857e+2,5.31278e+2,6.73482e+2,8.12590e+2,9.61242e+2])

sigma_vec1=np.array([1.89030e-40,8.07637e-41,3.01730e-41,1.23291e-41,5.26965e-42,2.25363e-42,1.15316e-42,5.15858e-43,3.45224e-43,2.30987e-43,1.03291e-43,4.83140e-44,2.16254e-44,9.67769e-45,4.33424e-45,1.94412e-45,9.54868e-46,6.70098e-46,5.38518e-46,6.49174e-46,7.81666e-46,9.84131e-46,1.29531e-45,1.63050e-45,2.45805e-45,3.23465e-45,4.25416e-45,5.59177e-45,6.71362e-45,7.70446e-45])

capture_vec_qm2=np.array([])
print('Sun, Spin-independent(qm2) interaction:')

for n,(mchi,sigma_p) in enumerate(zip(mchidd_vec1,sigma_vec1)):
    s1=timer()
    capture=WD.wimp_capture(WD.Sun, model_qm2, vmin, delta_eta, mchi, sigma_p=sigma_p, v_cut=18.5,targets_list=[WD.Al27])
    s2=timer()
    print('mchi=',mchi,'(',str((n+1)/len(mchidd_vec1)*100),'%) completed in ',s2-s1,' s')
    capture_vec_qm2=np.append(capture_vec_qm2,capture)
    
    
## In the work arXiv:1305.0912, qref=0.1 GeV is used (Table-1)
## so multiply wimpyc_capture with qref^2/qref^4=0.01
qref=0.1

model_qm4=WD.eft_hamiltonian('Spin-independent(qm4)',{(1,'qm4'): c_tau_SI_qm4})
WD.load_response_functions_capture(WD.Sun,model_qm4)

## Taken from arXiv:1305.0912
mchidd_vec2=np.array([4.09326e+0,4.30151e+0,4.47571e+0,4.79772e+0,5.35116e+0,6.02799e+0,6.79042e+0,7.80264e+0,8.96574e+0,1.03022e+1,1.19560e+1,1.45812e+1,1.81393e+1,2.32477e+1,3.16228e+1,4.21697e+1,6.02799e+1,8.11863e+1,1.12648e+2,1.57861e+2,2.16874e+2,2.92090e+2,3.89509e+2,5.09210e+2,6.72336e+2,8.19961e+2,9.80346e+2 ])

sigma_vec2=np.array([6.30014e-42,2.36944e-42,9.31498e-43,3.20625e-43,1.10445e-43,3.97832e-44,1.87153e-44,1.05234e-44,5.65966e-45,2.03944e-45,8.40017e-46,4.13785e-46,2.33023e-46,1.50053e-46,1.10551e-46,1.11167e-46,1.33743e-46,1.60718e-46,2.20840e-46,3.03512e-46,4.16972e-46,5.47706e-46,7.19291e-46,9.03175e-46,1.23985e-45,1.48707e-45,1.86401e-45])

capture_vec_qm4=np.array([])
print('Sun, Spin-independent(qm4) interaction:')

for n,(mchi,sigma_p) in enumerate(zip(mchidd_vec2,sigma_vec2)):
    s1=timer()
    capture=WD.wimp_capture(WD.Sun, model_qm4, vmin, delta_eta, mchi, sigma_p=sigma_p, v_cut=18.5,targets_list=[WD.Al27])
    s2=timer()
    print('mchi=',mchi,'(',str((n+1)/len(mchidd_vec2)*100),'%) completed in ',s2-s1,' s')

    capture_vec_qm4=np.append(capture_vec_qm4,capture)

pl.clf()
pl.plot(mchidd_vec1,capture_vec_qm2*qref**2,color='darkblue', linestyle='-',label='WimPyC, $q^{-2}$')
pl.plot(mchidd_vec2,capture_vec_qm4*qref**4,color='green', linestyle='-',label='WimPyC, $q^{-4}$')
pl.legend(frameon=False, fontsize=15)
pl.xscale('log')
pl.yscale('log')
pl.tick_params(axis='both', which='major', labelsize=16)  # Increase major tick labels
pl.tick_params(axis='both', which='minor', labelsize=12)
pl.tight_layout(pad=2.2)
pl.xlabel(r'$m_{\chi}$ [GeV]', fontsize=15)
pl.ylabel(r'$C$ [1/s]', fontsize=15)
pl.title(r'$C_1$, $^{27}$Al', fontsize=15)
tf=timer()
print('Total time taken:',tf-ti)
pl.show()



### Earth following 1609.08967
ti=timer()

# instantiates the effective Hamiltonian containing only O4/SD.    
WD.load_response_functions_capture(WD.Earth,model_SI)

# Halo function
vmin, delta_eta=WD.streamed_halo_function(v_esc_gal=533)

mchi_vec=np.logspace(1, 3, 20)

capture_vec=np.array([])
print('Earth, spin-dependent')

for n,mchi in enumerate(mchi_vec):
    s1=timer()
    capture=WD.wimp_capture(WD.Earth, model_SI, vmin=vmin, delta_eta=delta_eta, mchi=mchi, rho_loc=0.4, c0=1e-3/mv**2, c1=0,targets_list=[WD.P31])
    s2=timer()
    print('mchi=',mchi,'(',str((n+1)/20*100),'%) completed in ',s2-s1,' s')
    capture_vec=np.append(capture_vec,capture)

pl.clf()
pl.plot(mchi_vec,capture_vec,color='darkblue',label='WimPyC')
pl.legend(frameon=False, fontsize=15)
pl.xscale('log')
pl.yscale('log')
pl.xlim(mchi_vec[0],mchi_vec[-1])
pl.tick_params(axis='both', which='major', labelsize=16)  # Increase major tick labels
pl.tick_params(axis='both', which='minor', labelsize=12)
pl.tight_layout(pad=2.2)
pl.xlabel(r'$m_\chi$ [GeV]', fontsize=15)
pl.ylabel(r'C [1/s]', fontsize=15)
pl.title(r'$C_' + str(1) + r'$, $^{31}$P', fontsize=15)
tf=timer()
print('Total time taken:',tf-ti)
pl.show()


###Jupiter following 2411.04435
ti=timer()

WD.load_response_functions_capture(WD.Jupiter,WD.SD)

#Halo function
vmin, delta_eta=WD.streamed_halo_function()

mchi_vec=np.logspace(np.log10(1e-1), np.log10(10), 20)
capture_vec=np.array([])
for n,mchi in enumerate(mchi_vec):
    s1=timer()
    capture=WD.wimp_capture(WD.Jupiter, WD.SD, vmin=vmin, delta_eta=delta_eta, mchi=mchi, rho_loc=0.4, sigma_p=1e-35, cn_over_cp=0)
    s2=timer()
    capture_vec=np.append(capture_vec,capture)
    print('mchi=',mchi,'(',str((n+1)/20*100),'%) completed in ',s2-s1,' s')

pl.clf()
pl.plot(mchi_vec,capture_vec,color='darkblue',label='WimPyC')
pl.legend(frameon=False, fontsize=15)
pl.xscale('log')
pl.yscale('log')
pl.xlim(mchi_vec[0],mchi_vec[-1])
pl.tick_params(axis='both', which='major', labelsize=16)  # Increase major tick labels                                                                                
pl.tick_params(axis='both', which='minor', labelsize=12)
pl.tight_layout(pad=2.2)
pl.xlabel(r'$m_\chi$ [GeV]', fontsize=15)
pl.ylabel(r'C [1/s]', fontsize=15)
pl.title(r'Jupiter, $\sigma=10^{-35}$ cm$^2$', fontsize=15)
tf=timer()
print('Total time taken:',tf-ti)
pl.show()

    
### Main-sequence Star follwing 2405.12267
ti=timer()


# instantiates the effective Hamiltonian containing only O4/SD. 
WD.load_response_functions_capture(WD.MS_Star,WD.SI)

# Halo function
vmin, delta_eta=WD.streamed_halo_function(v_rot_gal=np.array([0.,1000.,0.]), v_esc_gal=550)

sigma_p_vec=np.logspace(np.log10(1e-39), np.log10(1e-30), 20)

capture_vec=np.array([])
print('Main-sequence Star, spin-independent')

for n,sigma_p in enumerate(sigma_p_vec):
    s1=timer()
    capture=WD.wimp_capture(WD.MS_Star, WD.SI, vmin=vmin, delta_eta=delta_eta, mchi=0.1, rho_loc=2e12, sigma_p=sigma_p, cn_over_cp=1)
    s2=timer()
    capture_vec=np.append(capture_vec,capture)
    print('mchi=',mchi,'(',str((n+1)/len(sigma_p_vec)*100),'%) completed in ',s2-s1,' s')

c_geom=WD.wimp_capture_geom(WD.MS_Star,mchi=0.1,vmin=vmin,delta_eta=delta_eta,rho_loc=2e12)


pl.clf()
pl.plot(sigma_p_vec, capture_vec, color='darkblue',label='WimPyC')
pl.plot(sigma_p_vec, c_geom*np.ones_like(sigma_p_vec), color='red',label='WimPyC (Geo)')
pl.legend(frameon=False, fontsize=15)
pl.xscale('log')
pl.yscale('log')
pl.xlim(1e-38,1e-30)
pl.tick_params(axis='both', which='major', labelsize=16)  # Increase major tick labels
pl.tick_params(axis='both', which='minor', labelsize=12)
pl.tight_layout(pad=2.2)
pl.xlabel(r'$\sigma_{p\chi}$ [cm$^2$]', fontsize=15)
pl.ylabel(r'$C$ [1/s]', fontsize=15)
pl.title(r'MStar, $M=$' +str(WD.MS_Star.mass)+'$~M_{\odot}$', fontsize=15)
tf=timer()
print('Total time taken:',tf-ti)
pl.show()
