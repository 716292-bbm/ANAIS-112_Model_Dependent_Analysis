import numpy 

import ctypes
import numpy as np
import os



from Final.Librerias.GeneraRitmo.Libs.RAPIDD.rapidd.core import _crapidd, base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo
from Final.Librerias.GeneraRitmo.Libs.RAPIDD.rapidd.experiments import counts_bin_LZ

# madDM values, maybe we can get them from the card? 

mneutron = 0.94  
mproton = 0.938

### SHOULD USE VALUES FROM PARAM CARD ### 

m_u = 2.550000e-03
m_d = 5.040000e-03
m_s = 1.010000e-01
m_c = 1.270000e+00
m_b = 4.700000e+00
m_t = 1.720000e+02


def cn_diracDM_SI(alphas, mchi, q=0.0):


    result = np.zeros(15) # Initialize a list with 31 zeros

    result[0] = (0.063255 * alphas["b_even"] +
        0.0628148 * alphas["c_even"]  +
        0.036 * alphas["d_even"]  +
        2 * alphas["d_odd"]  +
        0.0413 * alphas["s_even"]  +
        0.015 * alphas["u_even"] +
        alphas["u_odd"] 
    )
    #return result* vev**2 /(4*mchi*mneutron)
    return result 


def cp_diracDM_SI(alphas, mchi,  q=0.0):
    result = np.zeros(15) # Initialize a list with 31 zeros

    result[0] = (0.0632493 * alphas["b_even"] +
        0.00507276 * alphas["b_odd"] +
        0.0628148 * alphas["c_even"] -
        0.0125906 * alphas["c_odd"] +
        0.032 * alphas["d_even"] +
        1.0063 * alphas["d_odd"] +
        0.00629531 * alphas["e_odd"] +
        0.00629531 * alphas["mu_odd"] +
        0.0413 * alphas["s_even"] +
        0.00629531 * alphas["s_odd"] +
        0.00629531 * alphas["tau_odd"] +
        0.017 * alphas["u_even"] +
        1.98741 * alphas["u_odd"]
    )
    #return result* vev**2 /(4*mchi*mneutron)
    return result 

def cn_diracDM_SD(alphas, mchi, q=0.0):
    
    result = np.zeros(15) # Initialize a list with 31 zeros
    qsq = q**2
    result[3] = (
        -3.588 * alphas["d_even"]
        + 0.0170931 * alphas["d_odd"]
        + 0.124 * alphas["s_even"]
        + 0.00014071 * alphas["s_odd"]
        + 1.504 * alphas["u_even"]
        - 0.00205568 * alphas["u_odd"]
    )

    result[5] = (
        (0.679556 + 2.58587 * qsq) * alphas["d_even"]
        + (-0.0124848 - 0.685272 * qsq) * alphas["s_even"]
        + (-0.667072 - 1.9006 * qsq) * alphas["u_even"]
        / ((0.0182187 + qsq) * (0.300153 + qsq))
    )

    #return result * vev**2 /(4*mchi*mneutron)
    return result 


def cp_diracDM_SD(alphas, mchi, q=0.0):
    result = np.zeros(15) # Initialize a list with 31 zeros
    qsq = q**2
    result[3] = (
        1.504 * alphas["d_even"]
        - 0.00439169 * alphas["d_odd"]
        + 0.124 * alphas["s_even"]
        + 0.00014071 * alphas["s_odd"]
        - 3.588 * alphas["u_even"]
        + 0.00800105 * alphas["u_odd"]
    )

    result[5] = (
        (-0.667072 - 1.9006 * qsq) * alphas["d_even"]
        + (-0.0124848 - 0.685272 * qsq) * alphas["s_even"]
        + (0.679556 + 2.58587 * qsq) * alphas["u_even"]
    ) / ((0.0182187 + qsq) * (0.300153 + qsq))

    #return result* vev**2 /(4*mchi*mneutron)
    return result

alphasSI = {
    "d_even": 0.0,
    "d_odd": 0.0,
    "u_even": 0.0,
    "u_odd": 0.0,
    "s_even": 0.0,
    "s_odd": 0.0,
    "c_even": 0.0,
    "c_odd": 0.0,
    "b_odd": 0.0, 
    "b_even": 0.0, 
    "e_even": 0.0,
    "e_odd" : 0.0, 
    "mu_odd": 0.0,
    "tau_odd":0.0
}

alphasSD = {
    "d_even": 0.0,
    "d_odd": 0.0,
    "u_even": 0.0,
    "u_odd": 0.0,
    "s_even": 0.0,
    "s_odd": 0.0,
    "c_even": 0.0,
    "c_odd": 0.0,
    "b_odd": 0.0, 
    "b_even": 0.0, 
    "e_even": 0.0,
    "e_odd" : 0.0, 
    "mu_odd": 0.0,
    "tau_odd":0.0
}

def set_alphas_rpd(alphasSI, alphasSD, mchi, q=0.0):
    #print(alphasSI)
    cp = cp_diracDM_SI(alphasSI, mchi, q) + cp_diracDM_SD(alphasSD, mchi, q)
    cn = cn_diracDM_SI(alphasSI, mchi, q) + cn_diracDM_SD(alphasSD, mchi, q)

    c0,c1=isofromneuc(cp,cn)
    
    #print(c0)
    #print(c1)
    for i in range(15):
        set_any_Ncoeff(c0[i], i+1, "p")
        set_any_Ncoeff(c1[i], i+1, "n")
    
    return


###################### 
# Codes for the on-the fly sigma calculation
######################

def gen_Pff_mdm(card):


    Pff = {# the quark -> proton form factor in SIe SIo SDe SDo order
            "alpha_d": np.array([mproton*card['SPd']/m_d, card['VPd'], card['AVPd'], card['SigPd'] ]),
            "alpha_u": np.array([mproton*card['SPu']/m_u, card['VPu'], card['AVPu'], card['SigPu'] ]) ,
            "alpha_s": np.array([mproton*card['SPs']/m_s, 0.0, card['AVPs'], card['SigPs'] ]),
            "alpha_c": np.array([(2/27)*(mproton/m_c)*card['SPg'], 0.0, 0.0, 0.0 ]),##
            "alpha_b": np.array([(2/27)*(mproton/m_b)*card['SPg'], 0.0, 0.0, 0.0 ]),
            "alpha_t": np.array([(2/27)*(mproton/m_t)*card['SPg'], 0.0, 0.0, 0.0 ]) #
    }

    return Pff


def gen_Nff_mdm(card):
    Nff = {# the quark -> neutron form factor in SIe SIo SDe SDo order
            "alpha_d": np.array([mneutron*card['SNd']/m_d, card['VNd'], card['AVNd'], card['SigNd'] ]),
            "alpha_u": np.array([mneutron*card['SNu']/m_u, card['VNu'], card['AVNu'], card['SigNu'] ]) ,
            "alpha_s": np.array([mneutron*card['SNs']/m_s, 0.0, card['AVNs'], card['SigNs'] ]),
            "alpha_c": np.array([(2/27)*(mneutron/m_c)*card['SNg'], 0.0, 0.0, 0.0 ]),##
            "alpha_b": np.array([(2/27)*(mneutron/m_b)*card['SNg'], 0.0, 0.0, 0.0 ]),
            "alpha_t": np.array([(2/27)*(mneutron/m_t)*card['SNg'], 0.0, 0.0, 0.0 ]) # 
    }

    return Nff


def c1_dirac_mdm(card, alphaq):

    PFF = gen_Pff_mdm(card)
    NFF = gen_Nff_mdm(card)

    alpha_data = alphaq

    c1p, c1n = 0.0, 0.0 
    for quark, values in alpha_data.items():
        #print(values[0:2] * pff[quark][0:2])
        c1p += (values[0:2] * PFF[quark][0:2]).sum()
        c1n += (values[0:2] * NFF[quark][0:2]).sum()

    return c1p, c1n

def c4_dirac_mdm(card, alphaq):
    PFF = gen_Pff_mdm(card)
    NFF = gen_Nff_mdm(card)
    alpha_data = alphaq


    c4p, c4n = 0.0, 0.0 
    for quark, values in alpha_data.items():
        #print(values[0:2] * pff[quark][0:2])
        c4p += (values[2:] * PFF[quark][2:]).sum()
        c4n += (values[2:] * NFF[quark][2:]).sum()

    return c4p, c4n




def sigmaSI_nucleon_mdm(card, alphaq, dmmass):
    c1p, c1n = c1_dirac_mdm(card, alphaq)
    #print('c1p = %.2e c1n: %.2e' % (c1p, c1n))
    
    mup = dmmass * mproton / (dmmass + mproton)
    mun = dmmass * mneutron / (dmmass + mneutron)


    sigp = 4 *mup**2 * c1p**2 / np.pi
    sign = 4 *mun**2 * c1n**2 / np.pi

    #print('SI sigma_dmp: %.2e GeV^{-2} : %.2e cm^2' % (sigp, sigp*0.0389e-26))
    #print('SI sigma_dmn: %.2e GeV^{-2} : %.2e cm^2' % (sign, sign*0.0389e-26))

    return sigp, sign

def sigmaSD_nucleon_mdm(card, alphaq, dmmass):
    c4p, c4n = c4_dirac_mdm(card, alphaq)
    #print('c4p = %.2e c4n: %.2e' % (c4p, c4n))
    mup = dmmass * mproton / (dmmass + mproton)
    mun = dmmass * mneutron / (dmmass + mneutron)


    sigp = 12 *mup**2 * c4p**2 / np.pi
    sign = 12 *mun**2 * c4n**2 / np.pi

    #print('SD sigma_dmp: %.2e GeV^{-2} : %.2e cm^2' % (sigp, sigp*0.0389e-26))
    #print('SD sigma_dmn: %.2e GeV^{-2} : %.2e cm^2' % (sign, sign*0.0389e-26))

    return sigp, sign






if __name__== '__main__':

    rhochi_p = 0.3
    mchi_p = 50
    alphasSI["u_odd"]=(1e-20)**(1/4)
    alphasSI["d_odd"]=(1e-20)**(1/4)

    set_alphas_rpd(alphasSI, alphasSD, mchi_p) 
    
    read_halo()
  
    result = counts_bin_LZ(rhochi_p, mchi_p, 3.0, 60.0)* (60/1000)*2 *0.9

    print(result)