import numpy as np
from scipy import optimize
from scipy.interpolate import interp1d
from scipy.stats import poisson
import csv


def get_simple_limit (coeff,mchi,counts) :
    if counts==0 : 
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        return coeff * np.sqrt(2.3/counts)



def poisson_likelihood_limit (coeff,mchi,counts,bkgrd,observed) :
    if counts==0:
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        
        N_90 = float(optimize.root(lambda mu: poisson.cdf(observed, mu) - 0.1, 5).x)-bkgrd
        return np.sqrt(coeff**2* (N_90/counts))

def test_statistic(c,dm,bkgrd,obs) :
    #print(dm, bkgrd, obs)
    return np.sum( -2 * obs * np.log(((dm*c)+bkgrd)/bkgrd) + 2 * (dm * c) )

def binned_poisson_likelihood_limit (coeff,mchi,counts,data, background, cl=2.706) :
    '''cl is the desired confidence level, 0.2707 is for Cl table 40.2 pdg statistics'''
    observed = data
    bkgrd = background
    
    if np.sum(counts)==0:
        return np.inf
    else:
        #crosssec = calc_xsec(mchi,coeff)
        
        c_90 = float(optimize.root(lambda c: test_statistic(c,counts,bkgrd,observed)-cl, 1).x)
        return np.sqrt(c_90 * coeff**2)