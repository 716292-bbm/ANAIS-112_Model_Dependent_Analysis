import ctypes

import numpy as np
rapidd = ctypes.CDLL('build/libRAPIDD.so')

from scipy import optimize
from scipy.interpolate import interp1d
from scipy.stats import poisson
import csv


########

rapidd.set_any_Ncoeff.argtypes = [ ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
rapidd.set_any_Ncoeff.restype = ctypes.c_void_p

rapidd.read_halo.argtypes = [ctypes.c_char_p]
rapidd.read_halo.restype = ctypes.c_void_p

rapidd.read_efficiency.argtypes = [ctypes.c_char_p] 
rapidd.read_efficiency.restype = ctypes.c_void_p

rapidd.counts_effres_bin_LZ.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, 
                                        ctypes.c_char_p, ctypes.c_char_p]

rapidd.counts_effres_bin_LZ.restype = ctypes.c_double

#######

rapidd_lib = b"../lib/"
rapidd_lib_str = "../lib/"

rapidd.read_halo(rapidd_lib + b"/halo_table.dat") ### read the halo function 
pi = 3.1415
print("here")

def reset_coefficients() :
    for i in range(16):
        rapidd.set_any_Ncoeff(0, i, b"p")
        rapidd.set_any_Ncoeff(0, i, b"n")


def isofromneuc(cp, cn):
    ''' Simply takes coeffs from p n basis to 0 1 basis '''
    return (cp+cn)/2, (cp-cn)/2

def calc_xsec_SI(mchi,coeff) :
    mN = 0.938
    mv = 246
    pi = 3.14159
    u  = (mchi*mN)/(mchi+mN)
    return (float(u**2 * coeff**2 * 0.0389e-26/ (pi * mv**4)))

def calc_xsec_SD(mchi, coeff): 
    mN = 0.938
    mv = 246
    pi = 3.14159
    u  = (mchi*mN)/(mchi+mN)
    return (float(3.0*u**2 * coeff**2 * 0.0389e-26 / (16* pi * mv**4 )))



def lindhard (Er) :
    k=0.166
    Z=54 # so only for xenon
    e=11.5*Er*(Z**(-7/3))
    g=3*e**0.15 + 0.7*e**0.6 + e
    return k*g/(1+k*g)


def test_statistic(c,dm,bkgrd,obs) :
    #print(dm, bkgrd, obs)
    return np.sum( -2 * obs * np.log(((dm*c)+bkgrd)/bkgrd) + 2 * (dm * c) )


@np.vectorize
def counts_LZ_efficiency_res(rho, mass, E1, E2, eff_file = rapidd_lib + b'/efficiency_tables/LZ_NR_2022.csv'):
    rapidd.read_efficiency(eff_file) ### read efficiency table 
    return rapidd.counts_effres_bin_LZ( rho, mass, E1, E2, b"None" , b"ISO" )



e_kevee22, data22 = np.loadtxt(rapidd_lib_str+ '/experiments/lz2022/lz2022-data.csv',delimiter=',',unpack=True)
e_kevee22, bkgrd22 = np.loadtxt(rapidd_lib_str+ '/experiments/lz2022/lz2022-bkgrd.csv',delimiter=',',unpack=True)


def binned_poisson_likelihood_limit (coeff,mchi,counts,data, background) :
    
    observed = data
    bkgrd = background
    
    if np.sum(counts)==0:
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        
        c_90 = float(optimize.root(lambda c: test_statistic(c,counts,bkgrd,observed)-2.706, 1).x)
        return np.sqrt(c_90 * coeff**2)

def lzlimit(mchi, op=1, fnfp=1., coeff=1e-3, e_kevee=e_kevee22, data=data22, bkgrd=bkgrd22) :
    reset_coefficients()
    
    E1 = np.linspace(0,17,52)[0:51]
    E2 = np.linspace(0,17,52)[1:]
    bins = (E1+E2)/2
    spacing = bins[1]-bins[0]
    
    ## want to convert binning to kevnr to calculate the dm
    x=np.linspace(0,200,1000)
    y=lindhard(x)

    f_lind = interp1d(x*y,y)
    E1_lind = E1/f_lind(E1)
    E2_lind = E2/f_lind(E2)
    E1_lind[0]=0
    
    ## set the coefficients for the correct operator
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    ## scale the bkgrds and data down so that there are 11 bkgrd events
    totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(counts_LZ_efficiency_res(0.3,mchi,E1_lind,E2_lind) * (60/1000)*2 *0.9) 
    lim=(binned_poisson_likelihood_limit(coeff, mchi, dm , totaldata*scaling, totalbkgrd*scaling) )
        
    return lim


def lzlimitproj(mchi, op=1, fnfp=1., coeff=1e-3, eff_file=rapidd_lib + b'/efficiency_tables/LZ_NR_2022.csv',
               e_kevee=e_kevee22, bkgrd=bkgrd22) :
    reset_coefficients()
    
    ## import bkgrds and sort out binning in kevee units
    #e_kevee, data = np.loadtxt('lznew/lz2022-data.csv',delimiter=',',unpack=True)
    #e_kevee, bkgrd = np.loadtxt('lznew/lz2022-bkgrd.csv',delimiter=',',unpack=True)
    E1 = np.linspace(0,17,52)[0:51]
    E2 = np.linspace(0,17,52)[1:]
    bins = (E1+E2)/2
    spacing = bins[1]-bins[0]
    
    ## want to convert binning to kevnr to calculate the dm
    x=np.linspace(0,200,1000)
    y=lindhard(x)

    f_lind = interp1d(x*y,y)
    E1_lind = E1/f_lind(E1)
    E2_lind = E2/f_lind(E2)
    E1_lind[0]=0
    
    ## set the coefficients for the correct operator
    cp = coeff; cn = fnfp*coeff
    c0,c1=isofromneuc(cp,cn)

    rapidd.set_any_Ncoeff(c0, op, b"p") # ci, i (operator number), p: proton and n:neutron  
    rapidd.set_any_Ncoeff(c1, op, b"n") # ci, i (operator number), p: proton and n:neutron

    ## scale the bkgrds and data down so that there are 11 bkgrd events
    #totaldata = data*spacing
    totalbkgrd = bkgrd*spacing
    scaling = 11/(np.sum(totalbkgrd))
    
    ## calculate the dm
    limarray = []
    #for mchi in masses: 
    dm=(counts_LZ_efficiency_res(0.3,mchi,E1_lind,E2_lind,eff_file) * 2 ) #(60/1000)*2 ) 
    lim=(binned_poisson_likelihood_limit(coeff, mchi, dm , totalbkgrd*scaling*1000/60, totalbkgrd*scaling*1000/60) )
        
    return lim



if __name__== '__main__':

    import matplotlib.pyplot as plt 

    mspace= np.geomspace(1e0,1e3,50)
    
    LZresultSI = np.zeros(np.shape(mspace))
    LZresultSD = np.zeros(np.shape(mspace))
    LZprojSI = np.zeros(np.shape(mspace))
    LZprojSD = np.zeros(np.shape(mspace))
    


    ## import bkgrds and sort out binning in kevee units
    

    for i in range(len(mspace)):
        LZresultSI[i] = calc_xsec_SI(mspace[i], lzlimit(mspace[i], op=1))
        LZresultSD[i] = calc_xsec_SD(mspace[i], lzlimit(mspace[i], op=4, coeff=1e1))
        
        
        LZprojSI[i] = calc_xsec_SI(mspace[i], lzlimitproj(mspace[i], op=1))
        LZprojSD[i] = calc_xsec_SD(mspace[i], lzlimitproj(mspace[i], op=4, coeff=1e1))
    


    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))  # 1 row, 2 columns

    ax1.loglog(mspace, LZresultSI, label='LZ-2022')
    
    ax1.loglog(mspace, LZprojSI, ls='--', label='LZ future')
    
    
    ax1.set_xlabel(r'$m_{\rm DM}\,\,\left[{\rm GeV}\right]$')
    ax1.set_ylabel(r'$\sigma_{N}^{\rm SI}\,\,\left[{\rm cm}^2\right]$')

    
    
    ax2.loglog(mspace, LZresultSD)



    ax2.loglog(mspace, LZprojSD, ls='--')
    
    ax2.set_xlabel(r'$m_{\rm DM}\,\,\left[{\rm GeV}\right]$')
    ax2.set_ylabel(r'$\sigma_{N}^{\rm SD}\,\,\left[{\rm cm}^2\right]$')
    
    ax1.legend()
    
    plt.savefig('limits.pdf', bbox_inches='tight', dpi=120)
    plt.show()
