#WimPyDD_PATH='WimPyDD/'
#EXPERIMENTS_PATH=WimPyDD_PATH+'Experiments'
#HALO_FUNCTIONS_PATH=WimPyDD_PATH+'Halo_functions'
#WIMP_CAPTURE_PATH=WimPyDD_PATH+'WimPyC'



'''
WimPyDD is a modular, object–oriented and customizable Python code
that calculates accurate predictions for the expected rates in WIMP
direct–detection experiments within the framework of
Galilean–invariant non–relativistic effective theory in virtually any
scenario, including inelastic scattering, an arbitrary WIMP spin and a
generic WIMP velocity distribution in the Galactic halo. Starting from
version 2.0.0 the additional module WimPyC extends WimPyDD to the
calculation of WIMP capture in celestial bodies. WimPyDD and WimPyC
exploit the factorization of the three main components that enter in
the calculation of direct detection signals: i) the Wilson
coefficients of the effective theory, that encode the dependence of
the signals on the theoretical parameters; ii) a response function
that, besides nuclear physics, for DD depends on the features of the
experimental detector (acceptance, energy resolution, response to
nuclear recoils) and for WIMP capture on the properties of the
celestial body; iii) a halo function that depends on the WIMP velocity
distribution. The three components are calculated and stored
separately for later interpolation, combining them together only in
the last step of the signal evaluation procedure. This makes the
phenomenological study of the direct detection scattering rate with
WimPyDD or of WIMP capture with WimPyC transparent and fast also when
the parameter space of the WIMP model has a large dimensionality.

WimPyDD can be downloaded from:

https://wimpydd.hepforge.org/
'''

try:
    import numpy as np
except:
    import pip
    pip.main(['install', 'numpy'])
    import numpy as np
try:
    import matplotlib.pyplot as pl
except:
    import pip
    pip.main(['install', 'matplotlib'])
    import matplotlib.pyplot as pl
import os
from inspect import getargspec, signature
import pickle
try:
    import scipy.special as sp
except:
    import pip
    pip.main(['install', 'scipy'])
    import scipy.special as sp
    
from functools import reduce
import sys, select
import scipy.optimize as op
import scipy.integrate as integrate
import subprocess
import inspect
from scipy.integrate import tplquad
import scipy.optimize as op
from scipy.stats import poisson,norm

from ..package import *



def get_response_functions_interf_capture(isotope_obj,n_coeff1,n_coeff2=None,coeff1_input=None,coeff2_input=None, eft_amplitude_squared=eft_amplitude_squared,eft_amplitude_squared1=eft_amplitude_squared1,eft_modifier=lambda q:1.,j_chi=0.5,outputfile=None,verbose=True,n_sampling=1000,force_recalculation=False,filename_suffix='',u_esc_max=2000,er_min_cut=1e-6,er_max_cut=None,er_sampling=None,original_tuple=None): 

    if n_coeff2 is None:
        n_coeff2=n_coeff1

    # puts to zero non-interfering cuplings combinations
    # (need to update vanishing_response_functions routine!)
    if get_short_coupling(convert_to_all_spins(n_coeff1))!=get_short_coupling(convert_to_all_spins(n_coeff2)) and not symmetric_interference(convert_to_all_spins(n_coeff1),convert_to_all_spins(n_coeff2)) and not non_symmetric_interference(convert_to_all_spins(n_coeff1),convert_to_all_spins(n_coeff2)):
        # all response functions are zero. fills in a tuple with the correct stucture containing all zeros

        response_functions_interf=vanishing_response_functions_capture(isotope_obj)
        if verbose:
            print('response functions for c_'+print_coupling(n_coeff1)+'-c_'+print_coupling(n_coeff2)+' are vanishing')
        return response_functions_interf

    if coeff1_input is not None:
        if coeff2_input is None:
            coeff2_input=coeff1_input

    # creates output folders (WimPyC, WimPyC/Response_functions and WimPyC/Response_functions_spin_***) if missing (and adds __init__.py files)

    if not os.path.isfile(os.getcwd()+'/'+WIMP_CAPTURE_PATH+'/__init__.py'):
        open(os.getcwd()+'/'+WIMP_CAPTURE_PATH+'/__init__.py','a').close()

    path_response_functions=WIMP_CAPTURE_PATH+'/Response_functions'
    if not os.path.isdir(path_response_functions):
       p = subprocess.call("mkdir -p "+path_response_functions, stdout=subprocess.PIPE, shell=True)
    if not os.path.isfile(os.getcwd()+'/'+path_response_functions+'/__init__.py'):
        open(os.getcwd()+'/'+path_response_functions+'/__init__.py','a+').close()

    path_spin=path_response_functions+'/spin_'+print_spin(j_chi)
    if not os.path.isdir(path_spin):
       p = subprocess.call("mkdir -p "+path_spin, stdout=subprocess.PIPE, shell=True)
    if not os.path.isfile(os.getcwd()+'/'+path_spin+'/__init__.py'):
        open(os.getcwd()+'/'+path_spin+'/__init__.py','a+').close()
    
    path=path_spin

    # check time stamp of custom nuclear response functions, if present                                                                                                                                               
    func_w_path=os.getcwd()+'/WimPyDD/Targets/Nuclear_response_functions'

    time_w_func=0
    try:
        time_w_func=max(time_w_func,os.stat(func_w_path+'/'+isotope_obj.symbol+'_func_w.py').st_mtime)
    except:
        pass

    # loop to define all needed differential response functions
    # both to save them as attributes and to run primitive_table
    response_function_interf_tau_tau_prime=np.array([])
    for tau in range(2):
        response_function_interf_tau=np.array([])
        for tau_prime in range(2):
            if coeff1_input is None:
                coeff1=set_c_coeff_interf(n_coeff1,tau)
            else:
                coeff1=coeff1_input

            if coeff2_input is None:
                coeff2=set_c_coeff_interf(n_coeff2,tau_prime)
            else:
                coeff2=coeff2_input

            r=np.array([])
            r=np.append(r,def_response_function_capture(coeff1,isotope_obj,eft=lambda coeff1,q,element,n_isotope,coeff2: eft_modifier(q)*eft_amplitude_squared(coeff1,q,element,n_isotope,coeff2,j_chi),coeff2=coeff2))
            r=np.append(r,def_response_function_capture(coeff1,isotope_obj,eft=lambda coeff1,q,element,n_isotope,coeff2: eft_modifier(q)*eft_amplitude_squared1(coeff1,q,element,n_isotope,coeff2,j_chi),coeff2=coeff2))
            r=np.append(r,def_response_function_capture(coeff1,isotope_obj,eft=lambda coeff1,q,element,n_isotope,coeff2: eft_modifier(q)*eft_amplitude_squared1_e(coeff1,q,element,n_isotope,coeff2,j_chi),coeff2=coeff2))
            r=np.append(r,def_response_function_capture(coeff1,isotope_obj,eft=lambda coeff1,q,element,n_isotope,coeff2: eft_modifier(q)*eft_amplitude_squared1_em1(coeff1,q,element,n_isotope,coeff2,j_chi),coeff2=coeff2))

            response_function_interf_tau=np.append(response_function_interf_tau,r).reshape(-1,*r.shape)
        response_function_interf_tau_tau_prime=np.append(response_function_interf_tau_tau_prime,response_function_interf_tau).reshape(-1,*response_function_interf_tau.shape)

    isotope_obj.response_function_interf=response_function_interf_tau_tau_prime.reshape(2,2,4)

    # gets the info. for capture the
    # response functions are attributes of the isotope_obj target, not of the
    # celestial body. So needs only to check nuclear response functions
    # (to be compared with info in pre-existing integrated response
    # functions files, if present)
    info_new=isotope_obj.func_w.__doc__

    vel_dep_string=np.array(['R0','R1','R1*E_R','R1/E_R'])

    ## standard name for output file                                                                                                                                                                                  
    if n_coeff1 is not None and outputfile is None:
        inputfile=isotope_obj.symbol+'_c_'+print_coupling(n_coeff1)+'_c_'+print_coupling(n_coeff2)+filename_suffix+'.npy'
        inputfile2=isotope_obj.symbol+'_c_'+print_coupling(n_coeff2)+'_c_'+print_coupling(n_coeff1)+filename_suffix+'.npy'
    else:
        if outputfile is None:
            inputfile='custom_response_functions_'+filename_suffix+'.npy'
        else:
            inputfile=outputfile
            
        inputfile2=inputfile

    flag_tau_tau_prime_flip=0


    # loads pre-existing response functions if present
    if os.path.isfile(path+'/'+inputfile) and not force_recalculation:
        if verbose:
            print('loading '+inputfile)
            print('type: print(np.load(\''+path+'/'+inputfile+'\',allow_pickle=True)[-1]) to get info about input used')

        with open(path+'/'+inputfile, 'rb') as file:
            response_functions_interf=pickle.load(file,encoding='latin1')

        info_old=response_functions_interf[-1]

        time_response_functions=os.stat(path+'/'+inputfile).st_mtime

        n_count=0
        if info_old !=info_new and verbose:
            print(20*'*')
            print('WARNING: some input changed:')
            for (string_old,string_new) in zip(info_old,info_new):
                if string_old!=string_new:
                    n_count+=1
                    print(10*'-'+str(n_count)+10*'-')
                    print('info from response function file '+inputfile+':')
                    print(string_old)
                    print('info from input isotope:')
                    print(string_new)
            #in the case of capture the only mismatch can be on the
            # nuclear response functions
            print(20*'*'+'Warning'+20*'*'+'\nThe documentation collected from the response functions of the isotope '+isotope_obj.symbol+' does not match with that contained in the file '+path+'/'+inputfile+'.\n If any modification in the nuclear response functions of '+isotope_obj.symbol+' is more recent than the file '+path+'/'+inputfile+' you need to recalculate it. If this is not the case consider to update the documentation in the nuclear response functions of '+isotope_obj.symbol+' to stop this warning.')
            print(20*'*')

        if time_w_func>time_response_functions:
            print('**********************************************')
            print('*******************WARNING********************')
            print('**********************************************')
            print('a custom nuclear form factor for one of the targets is more recent then '+inputfile)

            print('if custom nuclear response functions need not to be updated use update_time_stamp=True ')
            print('when calling load_response_functions_capture to avoid this warning in the future')
            
            if force_recalculation:
                reply=1
            else:
                reply=int(delayed_input(10,'Type 1 if you want to recalculate the response functions, 2 if you want to use the response functions stored in '+path+'/'+inputfile+', 0 otherwise (default answer is 2 after 10 seconds)',2))
            if reply==1:
                response_functions_interf=()

                for n_vel in range(4):
                    if verbose:
                        print('calculating '+vel_dep_string[n_vel])
                    response_functions_tau=()
                    for tau in range(2):
                        response_functions_tau_prime=()
                        for tau_prime in range(2):
                            if verbose:
                                print('calculating tau='+str(tau)+', tau_prime='+str(tau_prime))
                            if tau>tau_prime and n_coeff1==n_coeff2:
                                response_functions_tau_prime+=(response_functions_tau[0][1],)

                            else:

                                # core of calculation. calls primitive table with
                                # input an array of response functions r with shape
                                # r[n_element][tau,tau_prime,n_isotope]
                                ## N.B. loop over n_vel is external to
                                ## primitive_table, loops over
                                ## n_element,tau,tau_prime, n_isotope
                                ## are internal for capture each file
                                ## contains only
                                ## [n_vel][tau][tau_prime][0] and
                                ## [n_vel][tau][tau_prime][1], so the
                                ## input r shape should be only
                                ## r[tau,tau_prime] (loop over tau and
                                ## tau_prime internal to
                                ## primitive_table and loop over n_vel
                                ## external). element and isotope are
                                ## fixed to input arguments
                                
                                if isinstance(original_tuple,type(None)): 
                                    response_functions_tau_prime+=(primitive_table_interf_capture(isotope_obj,tau,tau_prime,isotope_obj.response_function_interf[:,:,n_vel],n_sampling=n_sampling,verbose=verbose,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut,er_sampling=er_sampling),)
                                else:
                                    response_functions_tau_prime+=(primitive_table_interf_capture(isotope_obj,tau,tau_prime,isotope_obj.response_function_interf[:,:,n_vel],n_sampling=n_sampling,verbose=verbose,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut,er_sampling=er_sampling,original_tuple=original_tuple[n_vel][tau][tau_prime]),) 
        
                        response_functions_tau+=(response_functions_tau_prime,)
                    response_functions_interf+=(response_functions_tau,)
                response_functions_interf+=(info_new,)

                # saves output table                                                                                                                                                                                  
                with open(path+'/'+inputfile, 'wb') as file:
                    pickle.dump(response_functions_interf, file)

            elif reply==2:
                pass

            elif reply==0:
                print('Response functions for '+isotope_obj.symbol+' have not been loaded and response_functions_capture attribute set to the empty tuple (). Also the info in the file '+inputfile+' has not been updated.')
                response_functions_interf=()

        # ignored because   flag_tau_tau_prime_flip always zero
        else: #strange here; Kang is curious who wrote this?
            if flag_tau_tau_prime_flip==1:
                flip_tau_tau_prime_in_response_functions_interf(exp,response_functions_interf)
    # if pre-existing response functions are not present asks to calculate them                                                                                                                                       
    if not os.path.isfile(path+'/'+inputfile) or force_recalculation:

        if force_recalculation:
            print(path+'/'+inputfile+' will be overwritten')
            reply=1
        else:
            if coeff1_input is None:
                reply=int(delayed_input(10,inputfile+' not available for '+isotope_obj.symbol+' in '+path+'. Type 1 if you want to calculate the response functions, 0 otherwise (default answer is 1 after 10 seconds)',1))
            else:
                reply=int(delayed_input(10,'Response functions not available for '+isotope_obj.symbol+' in '+path+'. Type 1 if you want to calculate the response functions, 0 otherwise (default answer is 1 after 10 seconds)',1))


        if reply==1:
            response_functions_interf=()

            for n_vel in range(4):
                if verbose:
                    print('calculating '+vel_dep_string[n_vel])
                response_functions_tau=()
                for tau in range(2):
                    response_functions_tau_prime=()
                    for tau_prime in range(2):

                        if verbose:
                            print('calculating tau='+str(tau)+', tau_prime='+str(tau_prime))

                        if tau>tau_prime and n_coeff1==n_coeff2:
                            response_functions_tau_prime+=(response_functions_tau[0][1],)
                        else:

                            # core of calculation. calls primitive table with
                            # input an array of response functions r with shape
                            # r[n_element][tau,tau_prime,n_isotope]
                            ## N.B. loop over n_vel is external to
                            ## primitive_table, loops over
                            ## n_element,tau,tau_prime, n_isotope
                            ## are internal for capture each file
                            ## contains only
			    ## [n_vel][tau][tau_prime][0] and
                            ## [n_vel][tau][tau_prime][1], so the
                            ## input r shape should be only
                            ## r[tau,tau_prime] (loop over tau and
                            ## tau_prime internal to
                            ## primitive_table and loop over n_vel
                            ## external). element and isotope are
                            ## fixed to input arguments

                            if isinstance(original_tuple,type(None)): 
                                response_functions_tau_prime+=(primitive_table_interf_capture(isotope_obj,tau,tau_prime,isotope_obj.response_function_interf[:,:,n_vel],n_sampling=n_sampling,verbose=verbose,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut,er_sampling=er_sampling),)  
                            else:
                                response_functions_tau_prime+=(primitive_table_interf_capture(isotope_obj,tau,tau_prime,isotope_obj.response_function_interf[:,:,n_vel],n_sampling=n_sampling,verbose=verbose,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut,er_sampling=er_sampling,original_tuple=original_tuple[n_vel][tau][tau_prime]),)  

                    response_functions_tau+=(response_functions_tau_prime,)
                response_functions_interf+=(response_functions_tau,)
            response_functions_interf+=(info_new,)

            # ignored
            if flag_tau_tau_prime_flip==1:
                flip_tau_tau_prime_in_response_functions_interf(exp,response_functions_interf)

            #saves output
            with open(path+'/'+inputfile, 'wb') as file:
                pickle.dump(response_functions_interf, file)

        else:
            print('Response functions for '+isotope_obj.symbol+' have not been calculated and response_functions_capture attribute set to the empty tuple ().')
            response_functions_interf=()

    return response_functions_interf


def vanishing_response_functions_capture(isotope):
    
    try:
        info=isotope.response_functions_capture[-1]
    except:
        info=()

    response_functions_capture=()

    for n_vel in range(4):
        response_functions_tau=()
        for tau in range(2):
            response_functions_tau_prime=()
            for tau_prime in range(2):
                response_functions_tau_prime+=(np.zeros(200).reshape(2,-1),)
            response_functions_tau+=(response_functions_tau_prime,)
        response_functions_capture+=(response_functions_tau,)

    response_functions_capture+=(info,)

    return response_functions_capture


class isotope(object):
    '''
    Initializes a single nuclear isotope extracting it from an 
    object belonging to the "element" class.
    Isotope objects are needed as input by the following routines:

    WD.load_response_functions_capture
    WD.plot_response_functions_capture

    If the "element" object is not passed as a parameter it 
    is searched among the set of elements that is pre-defined in WimPyDD:

    Al,Ar,C,Ca,F,Fe,Ge,H,He,I,Mg,N,Na,Ne,Ni,O,S,Si,W,Xe 

    that can also be listed typing list_elements().
    Requires as input a string with the symbol of the isotope, 
    that must match one of the entries of the "isotopes" attribute 
    of the element object.

    Example: Xe128=WD.isotope('128Xe') 

    successfully initializes the Xe128 isotope if the predefined 
    element "WD.Xe" exists, and if its WD.Xe.isotopes attribute 
    contains the string '128Xe':

    >>>'128Xe' in WD.Xe.isotopes
    True

    If the isotope belongs to an element that is not pre-defined 
    one needs first to initialize a new element object (see help(WD.element) 
    and Appendix D of I, Jeong et al., Computer Physics Communications 
    276(2022)108342 for instructions on how to do it) and then 
    pass it as an argument to the isotope class. 
    For instance:

    >>> Cl35=WD.isotope('35Cl')
    looking for 35Cl in list of pre-defined targets:
    35Cl not found, check isotope spelling. Returned empty isotope 
    object with success attribute equal to False

    >>>cl=WD.element('Cl')
    '35Cl' in cl.isotopes
    True
    >>>cl.name
    chlorine
    >>>Cl35=WD.isotope('35Cl',element=cl)
    35Cl found in chlorine
    --------------------------------------------
    Parameters:                                                                           
    
    ​- symbol(str) - A string that identifies an isotope. By convention 
    should start with an atomic mass number followed by the capitalized 
    symbol of the element. Must match one of the entries of the "isotopes" 
    attribute of the "element" object input (see below)

​    Example: H1=WD.isotope('1H') requires that the predefined element WD.H 
    contains the string '1H' in its isotopes attribute:  
    WD.H.isotopes -> array(['1H', '2H'])                       

    - element(class element)  - The object belonging to the element class 
    where the isotope is contained. If None, the code will search among the 
    pre-defined elements in WimPyDD (type WD.list_elements() to list them). 
    Default=None

    ​- verbose(bool)  - If True prints out details of the process that 
    searches for the isotope. Default=True
    --------------------------------------
    Attributes: 

​    - a(integer)	       - atomic mass numbers of the isotope
​    - element(element object)  - the element object from where the isotope has been extracted
    - itar(integer)	       - internal code to access the default Nuclear functions 
      coefficients (see element object)
    - mass(float)	           - mass of the isotope in GeV 
    - spin(float)	           - spin of the isotope 
    - symbol(str)		   - symbol of the isotope                                                  
    - z(int)		   - atomic number of the isotope    
    - response_functions_capture(dict)    - a dictionary containing 
    the tabulated response functions loaded by the load_response_functions_capture routine
    (see help(WD.load_response_functions_capture routine) for details).
    contains the corresponding tabulated 
    sigma_tilde used by WD.wimp_capture_rate for the calculation of 
    WIMP capture rate signals in celestial bodies. 
    For a definition of the integrated response functions  sigma_tilde see WimPyC paper(arXiv:2510.21185).

    '''

    def __init__(self, symbol,element=None,verbose=True): #strange indentaion here; Kang is curiour who wrote this?
        present=False
        self.success=True
        if element is None:
            if verbose:
                print('looking for '+symbol+' in list of pre-defined targets:')
            for e in element_list:
                if symbol in e.isotopes:
                    element=e
                    present=True
                    if verbose:
                        print(symbol+' found in predefined '+e.name+' element object')
        else:
            if symbol in element.isotopes:
                if verbose:
                    print(symbol+' found in '+element.name)
                    present=True

        if not present:
            self.success=False
            if verbose:
                print(symbol+' not found, check isotope spelling. Returned empty isotope object with success attribute equal to False')
            return

        self.element=element

        n_isotope=np.where(symbol==element.isotopes)[0][0]
        self.n_isotope=n_isotope
        self.mass=element.mass[n_isotope]
        self.spin=element.spin[n_isotope]
        self.a=element.a[n_isotope]
        self.z=element.z
        self.data_mod=element.data_mod[n_isotope]
        self.itar=element.itar[n_isotope]
        self.symbol=symbol
        self.func_w=element.func_w[n_isotope]
        self.response_functions_capture=response_functions({})
        self.diff_response_functions_capture=diff_response_functions({})

    def __str__(self):
        return "symbol "+self.symbol+", atomic number "+str(self.z)+", mass "+str(self.mass)+'\n'+'Nuclear form factor:'+self.func_w.__doc__


def primitive_table_interf_capture(isotope_obj,tau,tau_prime,response_function,n_sampling=100,verbose=True,u_esc_max=2000,er_min_cut=1e-6,er_max_cut=None,er_sampling=None,original_tuple=None): 

    path=WIMP_CAPTURE_PATH

    if er_sampling is None:
        er01=er_min_cut
        if er_max_cut==None:
            er02=2*isotope_obj.mass*(u_esc_max/300)**2
        else:
            er02=er_max_cut

        nstep=n_sampling

        if verbose:
            print('Er interval:',er01,er02)

        n1=int(nstep/2)
        er_vec1=np.logspace(np.log10(er01),np.log10(er02/20),n1)
        er_vec2=np.linspace(er02/20,er02,nstep-len(er_vec1)+2)[1:]
        er_vec=np.append(er_vec1,er_vec2)
        er_vec_out=er_vec[1:]

    else:
        er_vec=np.append(np.array([er_min_cut]),er_sampling)
        er_vec_out=er_sampling
        nstep=len(er_vec_out)

    if verbose:
        print('nstep=',nstep)
        print('calculating: ',isotope_obj.symbol)

    r=response_function[tau,tau_prime]
    r_int_vec=np.array([])
    n_show=max(int(nstep/10),1)

    if isinstance(original_tuple,type(None)):
        integral_vec=np.array([])
        for i,(er_min,er_max) in enumerate(zip(er_vec[:-1],er_vec[1:])):
            integral_vec=np.append(integral_vec,integrate.quad(r,er_min,er_max, epsrel=1.e-08, limit=10000)[0])
            if i%n_show == 0:
                if verbose:
                    print(str(int(10*i/nstep))+'0% evaluated')
        r_int_vec=np.cumsum(integral_vec)
        final=np.append(er_vec_out,r_int_vec).reshape(2,nstep)
        
        return final
    
    else:
        er_origin,r_int_origin=original_tuple
        
        if set(er_sampling).issubset(set(er_origin)):
            print('No new energy values to calculate. Existing response functions left unchanged.')
            return original_tuple
        
        r_diff1=np.append(r_int_origin[0],r_int_origin[1:]-r_int_origin[:-1])

        er12,nn=np.unique(np.append(er_origin,er_vec_out),return_index=True) #er_vec_out=er_sampling

        r_diff12=np.append(r_diff1,np.zeros(len(er_vec_out)))[nn]
        er12_sorted=np.sort(er12)
        r_diff12_sorted=r_diff12[np.argsort(er12)]
	
        id1=np.ones(len(er_origin))
        id2=2*np.ones(len(er_vec_out))
        id12=np.append(id1,id2)[nn]

        prod12=id12[1:]*id12[:-1]
        if er_origin[0]<=er_vec_out[0]:
            prod12=np.append(np.array([1]),prod12)
        else:
            prod12=np.append(np.array([2]),prod12)

        er_vec=np.append(np.array([er_min_cut]),er12_sorted)
        er_min=er_vec[:-1]
        er_max=er_vec[1:]

        mask12=prod12>1
        r_diff_prod=integrate.quad_vec(lambda t:r(er_min[mask12]+t*(er_max[mask12]-er_min[mask12]))*(er_max[mask12]-er_min[mask12]),0,1,epsrel=1.e-08,limit=10000)[0]
        r_diff12_sorted[mask12]=r_diff_prod

        r_int_out=np.cumsum(r_diff12_sorted)

        final=np.append(er12_sorted,r_int_out).reshape(2,len(er12_sorted))

        return final
    

def def_response_function_capture(coeff1,isotope_obj,eft=eft_amplitude_squared,coeff2=None):
    if coeff2 is None:
        coeff2=coeff1

    def response_function_capture(er):

        c_light=3.e5 # in km/sec
        hbarc2=0.389e-27 # in GeV^2 cm^2
        m_nucleus=isotope_obj.mass
        q=np.sqrt(2.*m_nucleus*er*1e-6)# in GeV
        
        element=isotope_obj.element
        n_isotope=isotope_obj.n_isotope

        eft_val=eft(coeff1,q,element,n_isotope,coeff2=coeff2)

        sigma_tilde=hbarc2*10**(-6)*2*m_nucleus/(4.*np.pi)*c_light**2*eft_val

        return sigma_tilde

    return response_function_capture


def load_response_functions_capture_isotope(isotope_obj,hamiltonian,j_chi=0.5,reset=False,verbose=True,update_time_stamp=False,n_sampling=1000,force_recalculation=False,u_esc_max=2000,er_min_cut=1e-6,er_max_cut=None,er_sampling=None,original_tuple=None,**args): 
    '''
Loads or calculates the WIMP-capture integrated response functions for the isotope object 
isotope_obj passed as an input. For the effective Hamiltonian hamiltonian and the WIMP spin j_chi 
updates the dictionary isotope_obj.response_functions with the entry {(hamiltonian,j_chi):r}, with
r a tuple containing the tabulated values.
  For each Wilson coefficients combination contained in hamiltonian.coeff_squared_list the 
tabulated response functions are first searched in the subfolder WimPyC/Response_functions. 
If missing they are calculated, saved in the subfolder and loaded to the dictionary. 

  In the updated dictionary r=isotope.response_functions_capture[hamiltonian, j_chi] is a tuple containing the tabulated response functions for all Wilson coefficient combinations. Each call to load_response_functions_capture_isotope adds a (hamiltonian,j_chi) entry to the dictionary unless reset=True.

N.B. Notice that the the response functions are an attribute of the 
nuclear isotopes and not of the celestial bodies. In other words, two celestial 
bodies containing the same isotope share the same response functions. 

    Structure of tuple r:
​        ------------------------------------------------------------------------                       
​        Setting:
                                                       
​        n_cicj=0,...,len(hamiltonian_obj.coeff_squared_list)-1 = entry of the 
        array hamiltonian_obj.coeff_squared_list containing couplings pairs                        
​        n_vel=0,1,2,3,4 corresponding to a=0,1,1E,1E^-1 (amplitude decomposition,       
        see help on the eft_amplitude_squared routine)

​        tau=0,1    (nuclear isospin)                                           
​        tau_prime=0,1 (nuclear isospin)

​        the tuple entry:                                                      

​        r[n_cicj][n_vel][tau][tau_prime][0]                                

​        contains a sampling of recoil energy values in keV
                                
​        and:                                                    

​        r[n_cicj][n_vel][tau][tau_prime][1]                             

​    contains the corresponding tabulated integrated response functions 
    sigma_tilde used by WD.wimp_capture_rate for the calculation of 
    WIMP capture rate signals in celestial bodies. 
    For a definition of the sigma_tilde response functions see WimPyC paper(arXiv:2510.21185)
------------------------------------------------------------------------------------
  The routine issues a warning and allows to recalculate existing response functions tables if they are older that the input files in the celestial body folder (densities.tab, mass. tab, radius.tab etc). Setting update_time_stamp=True updates the time stamps of all response functions tables to avoid such warning.
------------------------------------------------------------------------------------
Example:
	sun=celestial_body('Sun')
	o4_o6qm2=eft_hamiltonian('o4_o6',{4: func1, (6,'qm2'):func2})
	o4_o6qm2.coeff_squared_list ->
	[(4, 4), (4, (6,'qm2')), ((6,'qm2'), 4), ((6,'qm2'), (6,'qm2'))]
	j_chi=0.5

	Subfolder where the response functions are saved/read:
	WimPyDD/WimPyC/Response_functions/spin_1_2/
	
	The response functions tables are saved/read according to each target 
        element inside the celestial body
	by passing isotope_obj through load_response_functions_capture_isotope.
	i.e.
	1H_c_4_c_4.npy
	3He_c_4_c_6_qm2.npy
	18O_c_4_c_6_qm2.npy
	23Na_c_6_qm2_c_4.npy
	56Fe_c_6_qm2_c_6_qm2.npy
	...

-------------------------------------------------------------------------

  Input: 

​	- istope_obj	- object belonging to celestial_body class or isotope class.               

​	- hamiltonian	- object belonging to eft_hamiltonian class.

​	- j_chi(float)	- WIMP spin. Default: 0.5

​	- reset(bool)	- Empties the isotope.response_functions dictionary before adding the output. Default: False

​	- verbose(bool)	- If True prints out a list of the response function tables that are loaded or written. Default: True

​	- update_time_stamp(bool) - If true updates the time stamp of all response functions tables. Default: True

​	- n_sampling(int)	- The number of points of response functions sampling. Default: True

​	- force_recalculation(bool)	- If True delete saved files of response function before calculation. Default: False

​	- u_esc_max(float)	- maximum escape velocity inside the celestial body.

​	- er_min_cut(float)	- minimum recoil energy inside the celestial body.

​	- er_max_cut(float)	- maximum recoil energy inside the celestial body.

----------------------------------------------------------------------------                                    

  Output:                                                                      

  Adds a tuple containing the tabulated response function to the dictionary entry isotope.response_functions_capture[hamiltonian,j_chi]             

----------------------------------------------------------------------------
    '''

    if reset:
        isotope_obj.response_functions_capture=response_functions({})

    if isinstance(list(hamiltonian.wilson_coefficients.keys())[0],int) and j_chi>1:
        print('NO RESPONSE FUNCTION LOADED')
        print('Hamiltonian "'+hamiltonian.name+'" has NR couplings in the notation of PHYSICAL REVIEW C 89, 065501 (2014).\n For a WIMP spin higher that 1 use hamiltonian with couplings in all spins format (X, l,s).')
        return

    if isinstance(list(hamiltonian.wilson_coefficients.keys())[0],tuple):
        if isinstance(list(hamiltonian.wilson_coefficients.keys())[0][0],int) and j_chi>1:
            print('NO RESPONSE FUNCTION LOADED')
            print('Hamiltonian "'+hamiltonian.name+'" has NR couplings in the notation of PHYSICAL REVIEW C 89, 065501 (2014).\n For a WIMP spin higher that 1 use hamiltonian with couplings in all spins format (X, l,s).')
            return

    j_chi_out=print_spin(j_chi)
    response_functions_output=()
    diff_response_functions_output={}

    for n_coeff_num,(__c1__,__c2__) in enumerate(hamiltonian.coeff_squared_list):

        outputfile=isotope_obj.symbol+'_c_'+print_coupling(__c1__)+'_c_'+print_coupling(__c2__)+'.npy'
        if update_time_stamp:
            path=WIMP_CAPTURE_PATH+'/Response_functions/spin_'+print_spin(j_chi)+'/'
            if os.path.isfile(path+outputfile):
                p = subprocess.call("touch "+path+outputfile, stdout=subprocess.PIPE, shell=True)

        if isinstance(original_tuple,type(None)):
            response_functions_output+=(get_response_functions_interf_capture(isotope_obj=isotope_obj,n_coeff1=get_short_coupling(__c1__),n_coeff2=get_short_coupling(__c2__),eft_modifier=lambda q: hamiltonian.coeff_squared_q_dependence(q,__c1__=__c1__,__c2__=__c2__) ,j_chi=j_chi,verbose=verbose,outputfile=outputfile,n_sampling=n_sampling,force_recalculation=force_recalculation,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut,er_sampling=er_sampling),) 
        else:
            response_functions_output+=(get_response_functions_interf_capture(isotope_obj=isotope_obj,n_coeff1=get_short_coupling(__c1__),n_coeff2=get_short_coupling(__c2__),eft_modifier=lambda q: hamiltonian.coeff_squared_q_dependence(q,__c1__=__c1__,__c2__=__c2__) ,j_chi=j_chi,verbose=verbose,outputfile=outputfile,n_sampling=n_sampling,force_recalculation=force_recalculation,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut,er_sampling=er_sampling,original_tuple=original_tuple[n_coeff_num]),) 

        diff_response_functions_output[hamiltonian.name,j_chi,__c1__,__c2__]=isotope_obj.response_function_interf

    isotope_obj.diff_response_functions_capture.update(diff_response_functions(diff_response_functions_output))
    isotope_obj.response_functions_capture[hamiltonian.name,j_chi]=response_functions_output 



def plot_response_functions_capture(input_obj,model,j_chi=0.5,coeff_squared_list=None,style='',linewidth=1,tuple=None,n_vel_list=None,tau_list=None,tau_prime_list=None,rescaling_factor=1.,scatter_plot=False, scatter_plot_point_size=20, scatter_plot_color=None,er_min_cut=1e-6,er_max_cut=None,**args): 
    '''
input_obj can be a nuclear isotope or a celestial body.
Plots the response functions contained in isotope.response_functions[hamiltonian,j_chi]
----------------------------------------------------------------------------
Input:

  isotope or celestial body - object belonging either to the isotope class or to the celestial_body class

  model - object belonging to eft_hamiltonian class

  j_chi -  WIMP spin. Default: 0.5

  coeff_squared_list - A list of the Wilson coefficients combinations for which the response functions are plotted (among those in hamiltonian.coeff_squared_list). If None, all are plotted. Default: None

​                    Example: [(('Omega', 0, 0), ('Omega', 0, 0))]
	​                    [(4,4), (6,6)]

	        Example: if o4_o6.coeff_squared_list -> [(4, 4), (4, 6), (6, 4), (6, 6)],

		       plot_response_functions_capture(Sun.targets[0],o4_o6,j_chi=0.5)
		       total of 64 response functions
		       16 non-vanishing response functions plotted

	        plot_response_functions(Sun.targets[0],o4_o6,n_vel_list=[0],tau_list=[0],tau_prime_list=[1], coeff_squared_list=[[6,6]])
	        plots a single response function for c_6,c_6, n_vel=0 (velocity-independent term in squared amplitude), tau=0, tau_prime=1
----------------------------------------------------------------------------------
  Output:

  Figure showing response functions for defined experiment exp as a function of recoil energy

----------------------------------------------------------------------------------- 
    '''

    if type(input_obj).__name__=='isotope':
        isotope_list=[input_obj]
    elif type(input_obj).__name__=='celestial_body':
        isotope_list=input_obj.targets

    for isotope_obj in isotope_list:
        plot_response_functions_capture_isotope(isotope_obj,model,j_chi=j_chi,coeff_squared_list=coeff_squared_list,style=style,linewidth=linewidth,tuple=tuple,n_vel_list=n_vel_list,tau_list=tau_list,tau_prime_list=tau_prime_list,rescaling_factor=rescaling_factor,scatter_plot=scatter_plot, scatter_plot_point_size=scatter_plot_point_size, scatter_plot_color=scatter_plot_color,er_min_cut=er_min_cut,er_max_cut=er_max_cut,**args) 


def plot_response_functions_capture_isotope(isotope,model,j_chi=0.5,coeff_squared_list=None,style='',linewidth=1,tuple=None,n_vel_list=None,tau_list=None,tau_prime_list=None,rescaling_factor=1.,scatter_plot=False, scatter_plot_point_size=20, scatter_plot_color=None,er_min_cut=1e-6,er_max_cut=None,**args): 
    '''
Plots the response functions contained in isotope.response_functions[hamiltonian,j_chi]
----------------------------------------------------------------------------
Input:

    isotope: object belonging to isotope class

    model: object belonging to eft_hamiltonian class

    j_chi: WIMP spin

    coeff_squared_list (default: None): A list of the Wilson coefficients combinations for which the response functions are plotted (among those in hamiltonian.coeff_squared_list).
                       By default all are plotted.

                       Example: [(('Omega', 0, 0), ('Omega', 0, 0))]
                                [(4,4), (6,6)]

                       Example: if o4_o6.coeff_squared_list -> [(4, 4), (4, 6), (6, 4), (6, 6)],

                               plot_response_functions_capture(Sun.targets[0],o4_o6,j_chi=0.5)
                               total of 64 response functions
                               16 non-vanishing response functions plotted

                       plot_response_functions(Sun.targets[0],o4_o6,n_vel_list=[0],tau_list=[0],tau_prime_list=[1], coeff_squared_list=[[6,6]]
                       plots a single response function for c_6,c_6, n_vel=0 (velocity-independent term in squared amplitude), tau=0, tau_prime=1
----------------------------------------------------------------------------------
    Output:

    Figure showing response functions for defined experiment exp as a function of recoil energy
-----------------------------------------------------------------------------------
    '''

    if coeff_squared_list is None:
        coeff_squared_list=model.coeff_squared_list

    if tuple is None:
        if (model.name,j_chi) in isotope.response_functions_capture.keys(): 
            sigma_vec=isotope.response_functions_capture[model.name,j_chi] 
        else:
            print('response functions not yet loaded for '+isotope.symbol)
            answer=input('do you want to load them? (y/n)\n')
            if answer=='y':
                load_response_functions_capture(isotope,model,j_chi,er_min_cut=er_min_cut,er_max_cut=er_max_cut)
                sigma_vec=isotope.response_functions_capture[model.name,j_chi] 
            else:
                print('nothing to plot')
                return
    else:
        sigma_vec=tuple

    if n_vel_list is None:
        n_vel_list=range(4)
    if tau_list is None:
        tau_list=range(2)
    if tau_prime_list is None:
        tau_prime_list=range(2)

    ntot=0
    n_non_vanishing=0
    for __c1__,__c2__ in coeff_squared_list:
        n=model.coeff_squared_list.index((__c1__,__c2__))
        for n_vel in n_vel_list:
            for tau in tau_list:
                for tau_prime in tau_prime_list:
                    er=sigma_vec[n][n_vel][tau][tau_prime][0]
                    sigma=sigma_vec[n][n_vel][tau][tau_prime][1]
                    ntot+=1

                    if max(abs(sigma))!=0:
                        n_non_vanishing+=1

                        if scatter_plot:
                            pl.scatter(er,rescaling_factor*sigma,s=scatter_plot_point_size,color=scatter_plot_color)
                        else:
                            pl.plot(er,rescaling_factor*sigma,style,linewidth=linewidth,**args)

    print('total of '+str(ntot)+' response functions')
    print(n_non_vanishing,' non-vanishing response functions plotted')

    
def conditional_input(argument,argument_name,path,filename,verbose=True):

    file=os.path.isfile(path+filename)

    if argument is None and not file:
        argument=input('Please input the celestial body '+argument_name+' in solar units.\n')
        f=open(path+filename,'w')
        help='#value input by command line\n'
        f.write(help)
        f.write(str(argument))
        if verbose:
            print('saved '+argument_name+' '+str(argument)+' in solar units in '+path+filename)
        f.close()
    elif argument is not None and not file:
        f=open(path+filename,'w')
        help='#value input by command line\n'
        f.write(help)
        f.write(str(argument))
        if verbose:
            print('saved '+argument_name+' '+str(argument)+' in solar units in '+path+filename)
        f.close()
    elif argument is None and file:
        f=open(path+filename,'r')
        help=''
        for line in f:
            try:
                argument=float(line)
            except:
                help=line
        f.close()
        if verbose:
            print('Loaded '+argument_name+' '+str(argument)+' in solar units from '+path+filename)

    elif argument is not None and file:
        f=open(path+filename,'r')
        help=''
        for line in f:
            try:
                argument2=float(line)
            except:
                help=line
        f.close()

        if argument == argument2:
            print('Input '+argument_name+' '+str(argument)+' in solar units already present in '+path+filename)
        else:
            print(argument_name+' '+str(argument)+' passed as argument not matching with value '+str(argument2)+' contained in '+path+filename+'. Which one do you want to use?\n')
            print('1) '+str(argument)+' (input value)\n2) '+str(argument2)+' (value in '+filename+')')
            answer=input('Please input 1 or 2:  ')
            if answer=='1':
                f=open(path+filename,'w')
                help='#value input by command line\n'
                f.write(help)
                f.write(str(argument))
                print('Input '+argument_name+' '+str(argument)+' in solar units overwritten in '+path+filename+')')
            elif answer=='2':
                print('Input '+argument_name+' '+str(argument)+' ignored, value '+str(argument2)+' in solar units from '+filename+' used')
                argument=argument2
            else:
                print('Illegal answer. By default input '+argument_name+' '+str(argument)+' ignored, value '+str(argument2)+'from '+filename+' content used')

    return argument,help


class celestial_body(object):
    '''
Instantiates a celestial body object containing the astrophysical information 
required to calculate the WIMP capture rate: radius and mass of the system, target nuclei, densities. 
All the astrophysical information must be collected in a subfolder of the 
WimPyDD/Celestial_bodies directory matching the string name passed as the first argument. 
Such folder must contain the following user-provided files:

​	densities.tab

​	radius.tab

​	mass.tab

  Example:

​	sun=WD.celestial_body('Sun')

  Subfolder:

​	WimPyDD/WimPyC/Celestial_bodies/Sun

The mass and radius must be in solar units, so in the case of the Sun the two files  mass.tab and
radius.tab will contain 1.
The structure of the densities.tab file is an optional first line starting with "#" 
containing help information (returned by help(sun))
followed by headers and by a series of columns containing the radial dependence of nuclear densities.
The headers of the file must include 'r_vec', 'rho_tot' and at least one 
recognizable isotope name (see help(WD.isotope)). if 'v_esc' is missing it will be calculated using the 
radial density. In the r_vec, rho_tot and rho_i attributes the radius is normalized to 1 and 
the densities are normalized to a unit mass (i.e. integral(r^2*dr*rho_tot)=1) using the information
contained in radius.tab and mass.tab. So in densities.tab the units of the radius and 
of all densities are arbitrary (nuclear densities must only have the same units as the total 
density 'rho_tot'). The only dimensional quantity in the file is v_esc, which, if present, must be 
in km/s.

Example:
--------------------------------------------------
# Solar Model AGSS09ph, A. Serenelli et al, Astrophys. J. Lett. 705 (2009) L123–L127, [0909.2668].
r_vec   rho_tot     1H                    4He               .....     v_esc
0.0015  150.5      74.63614626654312     130.65634553820425 .....  1339.5764710930496  
0.002   150.5      67.19046650038389     117.57028777936976 .....  1339.6029805466728  
.....
0.9815  0.001162  0.0009100800407837143  0.00028441817946676305 .....  617.4823945410445  
------------------------------------------------------------------------------------------------- 
​    Parameters:

​	- name(str) 	- name of celestial body. Corresponds to the name of 
                        the subfolder in the Celestial_bodies folder containing 
                        the files with the astrophysical input.

​	- verbose(bool)    - If True print out the details of the loading procedure from the 
                           input directory. Default: False

----------------------------------------------------------------------------------------------------
  Attributes:

​	- name(str)	- celestial body name.

​	- mass(float)	- celestial body mass in solar units.

​	- r_vec(array)	- radius interval array for densities of isotopes (normalized to 1). 

​	- radius(float)	- celestial body radius in solar units.

​	- rho_i(dict)	- dictionary of densities for each isotope normalized to the 
                          isotope total mass fraction.
                          Example: sun.rho_i['4He'] contains the radial density profile of Helium 4
                          for each value of the radius contained in sun.r_vec
                          integral(sun.r_vec, 4*np.pi*sun.rho_tot*sun.r_vec**2)->0.26
                          

​	- rho_tot(array)	- density profile of the total radial density of the celestial_body 
                                  for each value of the radius contained in sun.r_vec
                                  normalized to a unit total mass:
                                  integral(sun.r_vec, 4*np.pi*sun.rho_tot*sun.r_vec**2)->1

​	- target_names(array)	- an array containing the names of the target isotopes inside the 
                                   celestial body (loaded from the headers of the densities.tab file).

​	- targets(array)	- an array of objects belonging to the isotope class corresponding 
                                  to the targets inside celestial body.

​	- v_esc(array)	-         an array of escape velocities for each value of the radius 
                                  contained in sun.r_vec. In km/sec. 
    '''

    def __init__(self, name,mass=None, radius=None, T_c=None, verbose=False): 
        self.name=name
        path=WIMP_CAPTURE_PATH+'/Celestial_bodies/'+self.name+'/'

        # get first the mass and the radius                                                                                                                                                                           
        if not os.path.isdir(path):
            p = subprocess.call("mkdir -p "+path, stdout=subprocess.PIPE, shell=True)
        if not os.path.isfile(os.getcwd()+'/'+WIMP_CAPTURE_PATH+'/Celestial_bodies'+'/__init__.py'):
            open(os.getcwd()+'/'+WIMP_CAPTURE_PATH+'/Celestial_bodies'+'/__init__.py','x').close()
        if not os.path.isfile(os.getcwd()+'/'+WIMP_CAPTURE_PATH+'/Celestial_bodies/'+self.name+'/__init__.py'):
            open(os.getcwd()+'/'+path+'/__init__.py','x').close()

        #If densities.tab is missing sets a template celestial body consisting of a hydrogen sphere of constant density
        #updated celestial_body class with use of delayed_input when initialization forlder is missing
        if not os.path.isfile(os.getcwd()+'/'+path+'/densities.tab'):
            if verbose:
                print('densities.tab is missing in '+path)
                print('Sphere of constant density with single nuclear target will be created')
            target=delayed_input(10,'choose target (default: 1H)','1H')
            f=open(path+'/densities.tab','w')
            n_points=10
            f.write('#template celestial body, hydrogen sphere of constant density\n')
            f.write('r_vec     rho_tot     1H\n')
            [f.write(str(n/n_points)+'         1        1\n') for n in range(n_points+1)]
            f.close()


        mass,info=conditional_input(mass,'mass',path,'mass.tab',verbose=verbose)
        self.mass=documented_float(mass)
        self.mass.info=info

        radius,info=conditional_input(radius,'radius',path,'radius.tab',verbose=verbose)
        self.radius=documented_float(radius)
        self.radius.info=info

        T_c,info=conditional_input(T_c,'core_temperature',path,'core_temperature.tab',verbose=verbose)
        self.T_c=documented_float(T_c)
        self.T_c.info=info

        ### the radius and densities are loaded in arbitrary units and
        ### returned normalized to 1
        if os.path.isfile(WIMP_CAPTURE_PATH+'/Celestial_bodies/'+name+'/densities.tab'):
            help_line,r_vec,rho_tot,v_esc,rho_i=get_celestial_body(WIMP_CAPTURE_PATH+'/Celestial_bodies/'+name+'/densities.tab')

            self.densities_help_line=help_line

            self.r_vec=r_vec/max(r_vec) # radius normalized to 1
            
            rho_average=4*np.pi*integrate.simps(rho_tot*self.r_vec**2,self.r_vec) # density normalization 
            self.rho_average=rho_average 
            self.rho_tot=rho_tot/rho_average # volume integral of self.rho_tot normalized to 1
            
            if v_esc is None:
                print('WARNING: v_esc label not found in densities.tab file. Check spelling if appropriate. Calculating v_esc from total density')
                self.v_esc=get_celestial_body_v_esc(self.mass,self.radius,self.r_vec,self.rho_tot)
                default='no'
                string='Want to add the escape velocity to '+WIMP_CAPTURE_PATH+'/Celestial_bodies/'+name+'/densities.tab'+'? (yes,no, default=no)'
                reply=delayed_input(10,string,default)
                if reply=='yes':
                    add_v_esc_to_densities(self)
            else:
                self.v_esc=v_esc# in km/sec

            target_names=np.array(list(rho_i.keys()))
            targets=[]

            missing_targets=[]
            present_targets=[]
            
            self.rho_i={}
            self.mass_fractions={}
            self.nt={}
            for target_name in target_names:
                target=isotope(target_name,verbose=verbose)
                if target.success:
                    self.rho_i[target_name]=rho_i[target_name]/rho_average
                    self.mass_fractions[target_name]=integrate.simps(4*np.pi*self.r_vec**2*self.rho_i[target_name],self.r_vec)
                    targets.append(target)
                    present_targets.append(target_name)
                    solar_mass_in_GeV=1.116e57
                    self.nt[target_name]=self.mass_fractions[target_name]*self.mass*solar_mass_in_GeV/target.mass
                else:
                    missing_targets.append(target_name)

            self.targets=targets
            self.target_names=[t.symbol for t in targets]
            self.target_indices={t.symbol:n for n,t in enumerate(self.targets)}

            if verbose:
                print(20*'=')
                print('WARNING')
                print('Only some of the isotopes contained in')
                print(WIMP_CAPTURE_PATH+'/'+self.name+'/densities.tab')
                print('were loaded:')
                print(present_targets)
                print('the following targets could not be initialized and were ignored:')
                print(missing_targets)
                print(20*'=')

        else:
            print('Densities not provided. You should provide '+path+'densities.tab.')

    def __str__(self):
        path=WIMP_CAPTURE_PATH+'/Celestial_bodies/'+self.name+'/'
        return str(self.name)+", radius "+str(self.radius)+" (solar units), mass "+str(self.mass)+" (solar units), densities of following targets:\n"+str(self.target_names)+"\ncontained in\n"+path+'densities.tab.\n'+self.densities_help_line

    
def get_mass(file,verbose=True):
    if os.path.isfile(file+'.tab'):
        if verbose:
            print('loading '+file+'.tab')

        f = open(file+'.tab', 'r')
        output=()
        info='No help provided'

        for line in f:
            try:
                exec('get_mass.t='+line)
                if isinstance(get_mass.t, int) or isinstance(get_mass.t,float) or isinstance(get_mass.t,list) or isinstance(get_mass.t,np.ndarray):
                    output+=get_mass.t,
                else:
                    info=str(line)
            except:
                info=str(line)

        output+=info,

        return output
    else:
        if verbose:
            print(file+'.tab missing, using 1. as default. ')
        return (1,'No input file provided, using 1. as default')

def get_radius(file,verbose=True):
    if os.path.isfile(file+'.tab'):
        if verbose:
            print('loading '+file+'.tab')

        f = open(file+'.tab', 'r')
        output=()
        info='No help provided'

        for line in f:
            try:
                exec('get_radius.t='+line)
                if isinstance(get_radius.t, int) or isinstance(get_radius.t,float) or isinstance(get_radius.t,list) or isinstance(get_radius.t,np.ndarray):
                    output+=get_radius.t,
                else:
                    info=str(line)
            except:
                info=str(line)

        output+=info,

        return output
    else:
        if verbose:
            print(file+'.tab missing, using 1. as default. ')
        return (1,'No input file provided, using 1. as default')


def get_celestial_body(filename):
    '''
    If the first line starts with "#" is kept as help, otherwise empty help is returned

    The headers of the file must include 'r_vec', 'rho_tot' and 'v_esc', and at least one recognizable isotope name
    '''
    f = open(filename, 'r')

    first_line=f.readline().strip()
    if first_line[0]=='#':
        help_line=first_line
        headers=f.readline().strip().split()
    else:
        help_line='No help provided for densities.tab content'
        headers=first_line.split()

    headers=np.array(headers)

    content=np.array([])
    for line in f:
        line = line.strip().split()
        columns = np.array([float(x) for x in line])
        nn=len(columns)
        content=np.append(content,columns).reshape(-1,nn)

    content=content.T

    #looks for required headers. If some is missing returns error

    required_headers=['r_vec', 'rho_tot']

    required_vectors=np.array([])
    for header in required_headers:
        if header in headers:
            nn=np.where(headers==header)
            required_vectors=np.append(required_vectors,content[nn][0]).reshape(-1,len(content[nn][0]))
            content=np.delete(content,nn,axis=0)
            headers=np.delete(headers,nn,axis=0)
        else:
            print('header '+header+" is missing. 'r_vec' and 'rho_tot' are required. No output returned")
            return

    r_vec=required_vectors[0]
    rho_tot=required_vectors[1]

    # renormalize r_vec in cm and rho_tot in grams/cm^3

    if 'v_esc' in headers:
        nn=np.where(headers=='v_esc')
        required_vectors=np.append(required_vectors,content[nn][0]).reshape(-1,len(content[nn][0]))
        content=np.delete(content,nn,axis=0)
        headers=np.delete(headers,nn,axis=0)
        v_esc=required_vectors[2]
    else:
        # the calling celestial body class will calculate it
        v_esc=None

    # builds the dictionary with the radial densities
    
    rho_i={k:v for k,v in zip(headers,content)}

    return help_line,r_vec,rho_tot,v_esc,rho_i


def load_response_functions_capture(input_obj,hamiltonian,j_chi=0.5,reset=False,verbose=True,update_time_stamp=False,n_sampling=1000,force_recalculation=False,u_esc_max=2000,er_min_cut=1e-6,er_max_cut=None,er_sampling=None,increase_sampling=False,**args): 
    '''
    For the effective Hamiltonian hamiltonian and the WIMP spin j_chi loads the differential response functions d_sigma_tilde_der and loads or calculates the tabulated response functions sigma_tilde (for their definition see WimPyC paper(arXiv:2510.21185)) used by wimp_dd_capture to calculate the WIMP capture rate in a celestial body.
    If input_obj belongs to the isotope class it updates the dictionaries input_obj.diff_response_functions and input_obj.response_functions for the single isotope.
    If input_obj belongs to the celestial_body class it updates the dictionaries of all the
    isotopes contained in input_obj.targets
-------------------------------------------------------------------------
  Input: 

​	- input_obj	- object belonging to celestial_body class or isotope class.

​	- hamiltonian	- object belonging to eft_hamiltonian class.

​	- j_chi(float)	- WIMP spin. Default: 0.5

​	- reset(bool)	- Empties the isotope.response_functions dictionary before adding the output. Default: False

​	- verbose(bool)	- If True prints out a list of the response function tables that are loaded or written. Default: True

​	- update_time_stamp(bool) - If true updates the time stamp of all response functions tables. Default: True

​	- n_sampling(int)	- The number of points of response functions sampling. Default: True

​	- force_recalculation(bool)	- If True delete saved files of response function before calculation. Default: False

​	- u_esc_max(float)	- maximum escape velocity inside the celestial body.

​	- er_min_cut(float)	- minimum recoil energy inside the celestial body.

​	- er_max_cut(float)	- maximum recoil energy inside the celestial body.

        - er_sampling            - used only if increase_sampling=True. A list of energy values in keV used to increase the sampling of the response functions. Default: None
    
        - increase_sampling       - if True the sampling of the response functions is increased using the energy values contained in the er_sampling parameter. Ignored if er_sampling parameter=None.
    
      
    
----------------------------------------------------------------------------
  Output:  
    
  if input_obj = isotope_obj:

             - add the entry {(hamiltonian_obj.name,j_chi,tau,tau_prime): func_array} to the dictionary isotope_obj.diff_response_functions_capture
             - adds the entry {(hamiltonian_obj.name,j_chi):r}  to the dictionary   isotope_obj.response_functions_capture
    
  if input_obj= celestial_body_obj repeats the same procedure for every target in celestial_body_obj.

    func_array: array of differential response functions with shape (2,2,4) corresponding to tau, tau_prime and n_vel with tau,tau_prime=0,1 the nuclear isospin and n_vel=0,1,2,3 corresponds to the amplitude decomposition 0,1,1E, 1E^{-1} (see help of eft_amplitude_squared).   

    
    Structure of tuple r:
​        ------------------------------------------------------------------------                       
​        Setting:
                                                       
​        n_cicj=0,...,len(hamiltonian_obj.coeff_squared_list)-1 = entry of the 
        array hamiltonian_obj.coeff_squared_list containing couplings pairs                        
​        n_vel=0,1,2,3 corresponding to a=0,1,1E,1E^-1 (amplitude decomposition,       
        see help on the eft_amplitude_squared routine)

​        tau=0,1    (nuclear isospin)                                           
​        tau_prime=0,1 (nuclear isospin)

​        the tuple entry:                                                      

​        r[n_cicj][n_vel][tau][tau_prime][0]                                

​        contains a sampling of recoil energy values in keV
                                
​        and:                                                    

​        r[n_cicj][n_vel][tau][tau_prime][1]                             

​    contains the corresponding tabulated integrated response functions 
    sigma_tilde used by WD.wimp_capture_rate for the calculation of 
    WIMP capture rate signals in celestial bodies. 
    For a definition of the sigma_tilde response functions see WimPyC paper(arXiv:2510.21185).

        N.B. Notice that the the response functions are an attribute of the 
        nuclear isotopes and not of the celestial bodies. In other words, two celestial 
        bodies containing the same isotope share the same response functions. 


    Example:
    having defined the isotope Al27 and the hamiltonian h:
    
    >>> print(Al27)
    symbol 27Al, atomic number 13, mass 25.137
    Nuclear form factor:Definition from Phys.Rev.C 89 (2014) 6, 065501 
    (e-Print: 1308.6288[hep-ph]) used for nuclear W functions as default

    >>> print(h)
    Hamiltonian name:h
    Hamiltonian:c_4()* O_4+c_5()* O_5+c_6()* O_6
    Squared amplitude contributions:
    O_4*O_4, O_4*O_6, O_5*O_4, O_5*O_5, O_6*O_4, O_6*O_6

    the command line:
    >>> WD.load_response_functions_capture(Al27, h, j_chi=0.5)
    calculates for a WIMP of spin=0.5 (default) the response functions 
    sigma_tilde for Al27 or loads them from 
    WimPyDD/WimPyC/Response_functions/spin_1_2 and stores them in the 
    Al27.response_functions_capture dictionary.

    The sigma_tilde corresponding to the c5*c4 combination of effective couplings, 
    to the amplitude component 0 and to the isospin combination [tau,tau_prime]=[0,1] is in: 
    
    >>> r=Al27.response_functions_capture[h.name,0.5]

    >>> n_cicj=2
    because:
    >>> h.coeff_squared_list[n_cicj=2] 
    (5, 4)

    n_vel=0 (with the correspondence 0->"0", 1->"1", 2->"1E", 3->"1E^-1", 
    see help(WD.eft_amplitude_squared))

    tau=0
    tau_prime=1

    
    energy=r[n_cicj][n_vel][tau][tau_prime][0]  
    sigma_tilde=r[n_cicj][n_vel][tau][tau_prime][1]  



----------------------------------------------------------------------------
    '''

    if type(input_obj).__name__=='isotope':
        isotope_list=[input_obj]
    elif type(input_obj).__name__=='celestial_body':
        isotope_list=input_obj.targets

    for isotope_obj in isotope_list:
        load_response_functions_capture_isotope(isotope_obj,hamiltonian,j_chi=j_chi,reset=reset,verbose=verbose,update_time_stamp=update_time_stamp,n_sampling=n_sampling,force_recalculation=force_recalculation,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut)

        if increase_sampling:
            if er_sampling is None:
                print('increase_sampling=True requires to pass er_sampling')
                print('er_sampling is missing, nothing done')
                return
            elif not er_sampling is None:
                if len(er_sampling)<1:
                    print('trying to increase sampling but er_sampling is empty, nothing done')
                    return

            if force_recalculation:
                print('when increase_sampling=True  and force_recalculation=True')
                print('cannot recalculate and increase sampling at the same time!')
                print('sampling increase is ignored')
                return
            tuple_hamiltonian1=isotope_obj.response_functions_capture[hamiltonian.name, j_chi]
            isotope_obj.response_functions_capture[hamiltonian.name+'(before sampling increase)', j_chi]=tuple_hamiltonian1
            load_response_functions_capture_isotope(isotope_obj,hamiltonian,j_chi=j_chi,reset=reset,verbose=verbose,update_time_stamp=update_time_stamp,n_sampling=n_sampling,force_recalculation=True,u_esc_max=u_esc_max,er_min_cut=er_min_cut,er_max_cut=er_max_cut,er_sampling=er_sampling,original_tuple=tuple_hamiltonian1)
            tuple_hamiltonian_12=isotope_obj.response_functions_capture[hamiltonian.name, j_chi]

            n1=len(tuple_hamiltonian1[0][0][0][0][0])
            n12=len(tuple_hamiltonian_12[0][0][0][0][0])
            print(20*'=')
            print('capture response functions sampling increase for '+isotope_obj.symbol+',')
            print(hamiltonian.name+' hamiltonian and WIMP spin='+str(j_chi)+':')
            print('original sampling:'+str(n1))
            print('final sampling :'+str(n12))
            print(20*'=')

            dump_response_functions_capture(isotope_obj,hamiltonian,j_chi,tuple_hamiltonian_12)



def wimp_capture(celestial_body_obj,hamiltonian,vmin,delta_eta,mchi,j_chi=0.5,delta=0,rho_loc=0.3,response_functions=None,targets_list=None,coeff_squared_list=None,reset_response_functions=False,print_contributions=False,tau_range=range(2),tau_prime_range=range(2),n_vel_list=range(4),sum_over_streams=True,verbose=False,v_cut=0,choose_zero=True,increase_sampling=False,response_functions_er_sampling=None,**args): 
    '''
Calculates the WIMP capture rate in a celestial body.

The routine load_response_functions_capture is called if isotope.response_functions 
does not contain the required response functions for hamiltonian and j_chi.
--------------------------------------------------------------------------------------
  Input:

​	- celestial_body_obj	- object belonging to celestial_body class

​	- hamiltonian		- object belonging to eft_hamiltonian class

​	- vmin			- an array containing a list of WIMP asymptotic speeds in km/s.

​	- delta_eta		- an array containing the contribution of each stream to the halo 
                                  function eta(v) in (km/sec)^-1
				  (the routine streamed_halo_function can be used to calculate it) 

​	- mchi			- WIMP mass in GeV

​	- j_chi			- WIMP spin. Default: 0.5

​	- delta			- Mass splitting for inelastic scattering in keV. Default: 0

​	- rho_loc		- local dark matter density in units of GeV/cm^3. Default: 0.3

​	- response_functions	- if passed, overrides the set of response functions. 
                                  Must be a tuple with the same format of 
                                  isotope.response_functions[hamiltonian,j_chi]. 
                                  Default: None

​	- targets_list		- a list of objects belonging to the isotope class to override 
                                  the target list contained in celestial_body_obj.targets. 
                                  If None the targets in celestial_body_obj.targets are used. 
                                  Default: None

        - coeff_squared_list    - a list of Wilson coefficients combinations used to calculate the capture rate (same as the attribute coeff_squared_list of the eft_hamiltonian class)
 
                                  Default: None

​	- reset_response_functions	- If True empties the dictionary isotope.response_functions 
                                          and reloads the response functions for hamiltonian and j_chi. 
                                          Default: False

​	- print_contributions	- If True, prints out the contributions to the capture rate from each 
                                  combination of the Wilson coefficients. After each call to 
                                  wimp_capture such contributions can be access via the dictionary:

		                  wimp_capture.dc_dvi[(target.symbol,[c2,c2])]

                                  with [c1,c2] any couplings combination contained in 
                                  hamiltonian.coeff_squared_list

                                  Default: False

​	- tau_range, tau_prime_range	- overrides the values of tau, tau_prime in the double sum 
                                          over nuclear isospins sum_{tau} sum_{tau_prime} R^{tau tau_prime}W^{tau tau_prime} that enters the calculation of the rate. It allows to calculate the contribution to the rate of specific isospin combinations. Default: range(2)

​	- n_vel_list				- a list of values 0,1,2,3 corresponding to a=0,1,1E,1E^-1 (amplitude decomposition.  see help on eft_amplitude_squared routine). Default: range(4)

​	- sum_over_streams	- if False does not sum over velocity streams and returns for each target isotopes the contributions to the rate from each stream, in an array with shape vmin.shape + r_vec.shape. Default: True

​	- verbose			- Passed to load_response_functions, if called. Default: False

​	- v_cut				- v_cut parameters that controls the maximal aphelium that the bound orbit of a WIMP after scattering must have in order to be captured. Introduced to take into account gravitational disturbances of other celestial bodies. For details, see 1308.5897, 1204.5120, 1305.0912, 2408.09658. See also help of er_capture routine.

                                         Default: 0

​	- choose_zero		- If True set dc_dvi as zero when E_max is smaller than E_min

                                 Default: True

        - increase_sampling       - if True the sampling of the response functions is increased using the energy values contained in the response_functions_er_sampling parameter.

                                 Default: False

        - response_functions_er_sampling     - list of energy values in keV used to increase the sampling of the response functions if increase_sampling is True. If None a list of energy values is obtained using the get_wimp_capture_response_functions_energy_sampling routine (see dedicated help).  

                                 Default: None                                   
---------------------------------------------------------------------------------------
  Output:

  	a value of capture rate
---------------------------------------------------------------------------------------
Example:
	sun=celestial_body('Sun')
	o1=eft_hamiltonian('o1',{1:lambda:[1,1]})
	load_response_functions_capture(sun,o1,j_chi=0.5)
	vmin,delta_eta0=streamed_halo_function() (time-independent component of halo function for Maxwellian with default
					      values of parameters - see help on streamed_halo_function)
	wimp_capture(sun,o1,vmin,delta_eta0,mchi=100)
	gives capture rate
---------------------------------------------------------------------------------------
    '''

    flatten=False
    
    if vmin[0]==0:
        vmin=vmin[1:]
        delta_eta=delta_eta[1:]

    input_args={}

    input_args.update({k:1 for k in {e for e in hamiltonian.global_arguments}-{t for t in list(args.keys())}
                       if k not in list(hamiltonian.global_default_args.keys()) and k!='mchi' and k!='delta'})

    input_args.update({k:v for k,v in hamiltonian.global_default_args.items() if k not in list(args.keys())})

    input_args.update(args)

    c_light=3e5
    solar_radius_in_cm=6.9634e10 #to convert from solar radii to cm                                                                                                                                                   
    solar_mass_in_gram=1.989e+33 #to convert from solar mass to gram                                                                                                                                                  
    celestial_body_radius=celestial_body_obj.radius*solar_radius_in_cm
    celestial_body_mass=celestial_body_obj.mass*solar_mass_in_gram

    r_vec=celestial_body_obj.r_vec

    v_esc=celestial_body_obj.v_esc


    GeV_per_gram=5.62e23
    lambda_i=delta_eta*vmin

    if targets_list is None:
        targets_list=celestial_body_obj.targets

    if coeff_squared_list is None:
        coeff_squared_list=hamiltonian.coeff_squared_list


    if tau_range is None:
        tau_range=range(2)
    if tau_prime_range is None:
        tau_prime_range=range(2)

    dc_dvi=np.zeros((len(r_vec),len(vmin)))

    dc_dvi_dict={} 

    for target in targets_list:
        if verbose:
            print(target.symbol)

        if (hamiltonian.name,j_chi) not in target.response_functions_capture.keys() or increase_sampling: 
            if verbose:
                print('Loading response functions for '+target.symbol+', '+hamiltonian.name+' and spin '+str(j_chi))
            if increase_sampling:
                if response_functions_er_sampling is None:
                    response_functions_er_sampling=get_wimp_capture_response_functions_energy_sampling(target, [celestial_body_obj], [mchi],[delta],vmin,**args) 
                
            load_response_functions_capture(target,hamiltonian,j_chi,verbose=verbose,increase_sampling=increase_sampling,er_sampling=response_functions_er_sampling,**args) 

        sigma_tilde_vec=target.response_functions_capture[hamiltonian.name,j_chi] 

        for n_coeff_squared,(__c1__,__c2__) in enumerate(coeff_squared_list):

            if verbose:
                print(__c1__,__c2__)

            if flatten:
                c1_c2_arguments=np.array([])
            else:
                c1_c2_arguments=np.unique(np.array([e for t in [hamiltonian.arguments[__c1__],hamiltonian.arguments[__c2__]] for e in t]))

            args_c1_c2={}
            for arg in c1_c2_arguments:
                if arg!='mchi' and arg!='delta':
                    args_c1_c2[arg]=input_args[arg]

            if 'mchi' in c1_c2_arguments: args_c1_c2['mchi']=mchi
            if 'delta' in c1_c2_arguments: args_c1_c2['delta']=delta

            dc_dvi_tau_tau_prime=np.zeros((len(r_vec),len(vmin)))
            for tau in tau_range:
                for tau_prime in tau_prime_range:

                    target_mass=target.mass
                    
                    # in the calculation densities are needed in grams/cm^3
                    # however in the celestial body the volume integral of the density
                    # of each target is normalized to the corresponding fraction mass

                    rho_target_grams_per_cm3=celestial_body_mass/celestial_body_radius**3*celestial_body_obj.rho_i[target.symbol]

                    # number density of target, in cm^-3
                    rho_target=rho_target_grams_per_cm3*GeV_per_gram/target_mass

                    mu=1./(1./mchi+1./target_mass)

                    vstar=get_vstar(target_mass, mchi, delta)


                    w=np.sqrt(v_esc.reshape(-1,1)**2+vmin**2)
                    w_larger_than_vstar=np.choose(w>vstar,[vstar,w])

                    e1_temp=er_min(target.mass, mchi, delta, w_larger_than_vstar)

                    ecap=er_capture(mchi,delta,vmin,v_cut=v_cut) 

                    e1=np.choose(e1_temp>ecap,[ecap,e1_temp])

                    e2=er_max(target.mass, mchi, delta, w_larger_than_vstar) 

                    n_v=len(vmin)
                    n_r=len(rho_target)
                    one=np.ones(n_r*n_v).reshape(n_r,n_v)

                    p_vdep=1e5*np.array([one/vmin,w_larger_than_vstar**2/c_light**2/vmin.reshape(1,-1)-delta*1e-6/mu/vmin,-one*target_mass/(2*mu**2*vmin),-one*(delta*1e-6)**2/(2*target_mass*vmin)]) 

                    for n_vel in n_vel_list:
                        er=sigma_tilde_vec[n_coeff_squared][n_vel][tau][tau_prime][0]
                        sigma_tilde=sigma_tilde_vec[n_coeff_squared][n_vel][tau][tau_prime][1]

                        diff_sigma_tilde=target.diff_response_functions_capture[hamiltonian.name, j_chi, __c1__, __c2__][tau,tau_prime,n_vel]

                        # Pre-allocate the sigma_matrix
                        sigma_matrix = np.zeros_like(e1)

                        # Mask where we should use the integral method
                        mask_int = e2 / e1 > 1 + 1e-2
                        mask_diff = ~mask_int

                        # Compute only where needed
                        if np.any(mask_int):
                            sigma_matrix[mask_int] = np.interp(e2[mask_int], er, sigma_tilde) - np.interp(e1[mask_int], er, sigma_tilde)

                        if np.any(mask_diff):
                            sigma_matrix[mask_diff] = 0.5 * (diff_sigma_tilde(e1[mask_diff]) + diff_sigma_tilde(e2[mask_diff])) * (e2[mask_diff] - e1[mask_diff])
                        

                        dc_dvi_target=rho_loc/mchi*rho_target.reshape(-1,1)*lambda_i.reshape(1,-1)*p_vdep[n_vel]*(sigma_matrix)*hamiltonian.coeff_squared(__c1__,__c2__,**args_c1_c2)[tau,tau_prime]
                        dc_dvi_target=np.choose(w>vstar,[0,dc_dvi_target])
                        coefficient=rho_loc/mchi*rho_target.reshape(-1,1)*lambda_i.reshape(1,-1)*p_vdep[n_vel]*hamiltonian.coeff_squared(__c1__,__c2__,**args_c1_c2)[tau,tau_prime]


                        if choose_zero:
                            dc_dvi_target=np.choose(e2>e1,[0,dc_dvi_target])
                            coefficient=np.choose(e2>e1,[0,coefficient]) 

                        dc_dvi+=dc_dvi_target
                        dc_dvi_tau_tau_prime+=dc_dvi_target


            dc_dvi_dict[target.symbol,(__c1__,__c2__)]=dc_dvi_tau_tau_prime 

    wimp_capture.dc_dvi=dc_dvi_dict 
    
    if sum_over_streams:
        capture=celestial_body_radius**3*np.sum(integrate.simps(4*np.pi*(dc_dvi.T*r_vec**2),r_vec))

    else:
        capture=celestial_body_radius**3*integrate.simps(4*np.pi*(dc_dvi.T*r_vec**2),r_vec)

    wimp_capture.capture=capture

    if print_contributions:
        print_wimp_capture_contributions(celestial_body_obj,dc_dvi) 

    return capture


def er_capture(mchi,delta,vmin,v_cut=0): 
    '''
Calculates the minimum energy that the WIMP needs to lose in order to be captured.
---------------------------------------------------------------------------------------

Input:

​	- mchi	- WIMP mass in GeV

​	- delta	- mass splitting in keV

​	- vmin	- incoming WIMP speed in km/s

​	- v_cut	-  in km/s. If different from zero counts as ‘‘captured’’ only dark matter particles 
                   which are kinematically constrained to orbits with maximum radius r0. 
                   Indicating with v_esc(r->infinity) the escape velocity from radius r and with 
                   v_esc(r->r0) the velocity to escape from radius r  to radius r0 
                   the velocity v_cut is given by the constant value
                   v_cut^2=v_esc(r->infinity)^2-v_esc(r->r0)^2. For instance, for r0=Sun-Jupiter
                   distance one has v_cut=18 km/s. Setting v_cut=18 km/s implies to 
                   count as captured only dark matter particles confined to orbits that lie 
                   within Jupiter’s orbit (see Eq.(8) of Kumar et al., PRD86 (2012) 073002, 
                   (arXiv:1204.5120)). 
                   Default: 0
  --------------------------------------------------------------------------------------
Output:
​	the value of minimum energy in keV that the WIMP needs to lose in order to be captured.
    '''
    return 1/2*mchi*((vmin/300)**2+(v_cut/300)**2)-delta 



def print_wimp_capture_contributions(celestial_body_obj,dc_dvi_contributions): 

    solar_radius_in_cm=6.9634e10
    celestial_body_radius=celestial_body_obj.radius*solar_radius_in_cm

    print('\n')

    tot={}
    capture_tot=0.
    for (target_name,(__c1__,__c2__)),dc_dvi in zip(dc_dvi_contributions.keys(),dc_dvi_contributions.values()):

        capture=celestial_body_radius**3*np.sum(integrate.simps(4*np.pi*(dc_dvi.T*celestial_body_obj.r_vec**2),celestial_body_obj.r_vec))

        tot[target_name,(__c1__,__c2__)]=capture
        
        capture_tot+=capture

        print(20*'=')
        print(__c1__,__c2__)
        print(target_name+':'+str(capture)+', '+str(capture/wimp_capture.capture*100)+'%')
        print(20*'=')
    print('total capture=',capture_tot)
    print_wimp_capture_contributions.totals=tot




def wimp_capture_matrix(celestial_body_obj,hamiltonian,vmin,delta_eta,mchi,delta=0,j_chi=0.5,targets_list=None,verbose=False,**args): 
    '''
    Outputs the matrix M that calculates the capture rate C for a couplings array c:
    C=np.dot(c,np.dot(M),c))
    for a given celestial body, hamiltonian, and velocity distribution.
    A dictionary with the mapping between the array components of c and the couplings of the hamiltonian is provided by the routine get_mapping.
    N.B. The output matrix M of wimp_dd_matrix is always in the isospin base. The M_pn in the proton-neutron base can be obtained using:
    U=WD.rotation_from_isospin_to_pn(hamiltonian)
    M_pn=np.dot(U,np.dot(M,U))
    ----------------------------------------------
    Input

    - celestial_body_obj - object belonging to the celestial_body class
    - hamiltonian - object belonging to eft_hamiltonian class
    - vmin - array containing a list of WIMP stream speed velocities in the lab frame in km/s.
    - delta_eta - array containing the contribution of each stream to the halo function eta(v) in (km/sec)^-1
                   (the routine streamed_halo_function can be used to calculate it)
    - mchi - WIMP mass in GeV
    - delta - Mass splitting for inelastic scattering in keV. Default: 0
    - j_chi - WIMP spin. Default: 0.5
    - targets_list -  a list of objects belonging to the isotope class to override the list of targets contained in celestial_body_obj.targets. 
                   If None the targets in celestial_body_obj.targets are used. 
                   Default: None

    - verbose - Passed to load_response_functions, if called. Default: False
    **args: any argument that can be passed to the wimp_capture routine
    ----------------------------------------------
    Output: array with shape (2*len(hamiltonian.couplings),2*len(hamiltonian.couplings))
    (the dimensionality corresponds to twice the number of couplings to include two values of the nuclear isospin)

    Example:

    >>>c_1_3=WD.eft_hamiltonian('c_1_3',{1: lambda: [1,1],3: lambda: [1,1]})
    >>>mapping=WD.get_mapping(c_1_3)

    >>>mapping[3,0]
    >>>2

    the index of the coupling c_3^0 is 2.

    >>>sun=WD.celestial_object('Sun')

    >>>m=WD.wimp_capture_matrix(sun,c_1_3,n_bin,vmin,delta_eta,mchi)
    '''

    solar_radius_in_cm=6.9634e10 #to convert from solar radii to cm
    celestial_body_radius=celestial_body_obj.radius*solar_radius_in_cm
    
    h=hamiltonian

    dim=2*len(hamiltonian.couplings)
    mapping=get_mapping(hamiltonian)
    matrix=np.zeros(dim*dim).reshape(dim,dim)

    if targets_list is None:
        targets_list=celestial_body_obj.targets

    for tau in range(2):
        for tau_prime in range(2):
            wimp_capture(celestial_body_obj, h, vmin, delta_eta, mchi,j_chi=j_chi, delta=delta,tau_range=[tau],tau_prime_range=[tau_prime],verbose=verbose,**args)
            for c1,c2 in hamiltonian.coeff_squared_list:
                m=mapping[c1,tau]
                m_prime=mapping[c2,tau_prime]
                capture=0
                for target in targets_list: 
                    symbol=target.symbol
                    dc_dvi=wimp_capture.dc_dvi[symbol,(c1,c2)]
                    capture+=celestial_body_radius**3*np.sum(integrate.simps(4*np.pi*(dc_dvi.T*celestial_body_obj.r_vec**2),celestial_body_obj.r_vec))

                matrix[m,m_prime]=capture
                if non_symmetric_interference(c1,c2):
                    matrix[m,m_prime]=matrix[m,m_prime]/2.

                matrix[m_prime,m]=matrix[m,m_prime]

    return (matrix+matrix.T)/2


def get_celestial_body_v_esc(mass,radius,r_vec,rho_tot):
    '''
    r_vec and rho_tot normalized to 1. mass and radius in solar units. output in km/sec
    '''

    if r_vec[0]==0:
        r_vec[0]=r_vec[-1]*1e-6


    solar_mass_in_grams=mass*1.989e33 
    solar_radius_in_cm=radius*6.9634e10
    r_vec_in_cm=solar_radius_in_cm*r_vec
    G=6.6743e-8 #cm^3*g^{-1}*s^{-2}
    rho_tot_grams_per_cm3=solar_mass_in_grams/solar_radius_in_cm**3*rho_tot


    func1=lambda r: r**2*np.interp(r,r_vec_in_cm,rho_tot_grams_per_cm3)
    func2=lambda r: r*np.interp(r,r_vec_in_cm,rho_tot_grams_per_cm3)

    r0=r_vec_in_cm[-1]

    v_esc_squared_1=np.array([])
    v_esc_squared_2=np.array([])

    for r in r_vec_in_cm:
        v_esc_squared_1=np.append(v_esc_squared_1,1e-10*8*np.pi*G/r*integrate.quadrature(func1, 0, r, maxiter=1000)[0])
        v_esc_squared_2=np.append(v_esc_squared_2,1e-10*8*np.pi*G*integrate.quadrature(func2, r, r0, maxiter=1000)[0])

    v_esc_squared=v_esc_squared_1+v_esc_squared_2

    v_esc=np.sqrt(v_esc_squared)

    return v_esc


def add_v_esc_to_densities(celestial_body_obj):

    file_name=WIMP_CAPTURE_PATH+'/Celestial_bodies/'+celestial_body_obj.name+'/densities.tab'

    #reads file
    f=open(file_name, 'r')
    input_content=list(f)
    f.close()

    #overwrites file adding v_esc in the last column
    f=open(file_name, 'w')

    nn=0
    for line in input_content:
        if line[0]=='#':
            f.write(line)
        elif 'r_vec' in line:
            f.write(line[:-1]+'  v_esc'+'\n')

        else:
            f.write(line[:-1]+'  '+str(celestial_body_obj.v_esc[nn])+'\n')
            nn+=1

    f.close()

    print('escape velocity added to '+file_name)


def d_wimp_capture_dvi_de(celestial_body_obj,hamiltonian,t,vmin,delta_eta,mchi,j_chi=0.5,delta=0,rho_loc=0.3,targets_list=None,coeff_squared_list=None,tau_range=range(2),tau_prime_range=range(2),n_vel_list=range(4),verbose=False,v_cut=0,choose_zero=True,**args):
    '''
    calculate dC/(dvdE) used by wimp_capture_accurate.
    ------------------------------------------------------
    Input:
    ------------------------------------------------------
    Output:
    '''
    if vmin[0]==0:
        vmin=vmin[1:]
        delta_eta=delta_eta[1:]

    input_args={}

    input_args.update({k:1 for k in {e for e in hamiltonian.global_arguments}-{t for t in list(args.keys())}
                       if k not in list(hamiltonian.global_default_args.keys()) and k!='mchi' and k!='delta'})

    input_args.update({k:v for k,v in hamiltonian.global_default_args.items() if k not in list(args.keys())})

    input_args.update(args)

    c_light=3e5
    solar_radius_in_cm=6.9634e10 #to convert from solar radii to cm

    solar_mass_in_gram=1.989e+33 #to convert from solar mass to gram

    celestial_body_radius=celestial_body_obj.radius*solar_radius_in_cm
    celestial_body_mass=celestial_body_obj.mass*solar_mass_in_gram

    r_vec=celestial_body_obj.r_vec

    v_esc=celestial_body_obj.v_esc

    GeV_per_gram=5.62e23
    lambda_i=delta_eta*vmin

    if targets_list is None:
        targets_list=celestial_body_obj.targets

    if coeff_squared_list is None:
        coeff_squared_list=hamiltonian.coeff_squared_list

    dc_dvi_de=np.zeros((len(r_vec),len(vmin)))

    if tau_range is None:
        tau_range=range(2)
    if tau_prime_range is None:
        tau_prime_range=range(2)

    for target in targets_list:
        dc_dvi_de_target=np.zeros((len(r_vec),len(vmin)))
        target_mass=target.mass

        if verbose:
            print(target.symbol)

        for n_coeff_squared,(__c1__,__c2__) in enumerate(coeff_squared_list):

            if verbose:
                print(__c1__,__c2__)

            c1_c2_arguments=np.unique(np.array([e for t in [hamiltonian.arguments[__c1__],hamiltonian.arguments[__c2__]] for e in t]))

            args_c1_c2={}
            for arg in c1_c2_arguments:
                if arg!='mchi' and arg!='delta':
                    args_c1_c2[arg]=input_args[arg]

            if 'mchi' in c1_c2_arguments: args_c1_c2['mchi']=mchi
            if 'delta' in c1_c2_arguments: args_c1_c2['delta']=delta

            dc_dvi_tau_tau_prime=np.zeros((len(r_vec),len(vmin)))
            for tau in tau_range:
                for tau_prime in tau_prime_range:

                    if verbose:
                        print('tau,tau_prime:',tau,tau_prime)

                    n_coeff1=get_short_coupling(__c1__)
                    n_coeff2=get_short_coupling(__c2__)
                    coeff1=set_c_coeff_interf(n_coeff1,tau)
                    coeff2=set_c_coeff_interf(n_coeff2,tau_prime)

                    eft_modifier=lambda q: hamiltonian.coeff_squared_q_dependence(q,__c1__=__c1__,__c2__=__c2__)

                    # in the calculation densities are needed in grams/cm^3
                    # however in the celestial body the volume integral of the density
                    # of each target is normalized to the corresponding fraction mass

                    rho_target_grams_per_cm3=celestial_body_mass/celestial_body_radius**3*celestial_body_obj.rho_i[target.symbol]

                    # number density of target, in cm^-3
                    rho_target=rho_target_grams_per_cm3*GeV_per_gram/target_mass

                    mu=1./(1./mchi+1./target_mass)

                    vstar=get_vstar(target_mass, mchi, delta)

                    w=np.sqrt(v_esc.reshape(-1,1)**2+vmin**2)
                    w_larger_than_vstar=np.choose(w>vstar,[vstar,w])

                    e1_temp=er_min(target.mass, mchi, delta, w_larger_than_vstar)

                    ecap=er_capture(mchi,delta,vmin,v_cut=v_cut) 

                    e1=np.choose(e1_temp>ecap,[ecap,e1_temp])

                    e2=er_max(target.mass, mchi, delta, w_larger_than_vstar) 

                    n_v=len(vmin)
                    n_r=len(rho_target)
                    one=np.ones(n_r*n_v).reshape(n_r,n_v)

                    p_vdep=1e5*np.array([one/vmin,w_larger_than_vstar**2/c_light**2/vmin.reshape(1,-1)-delta*1e-6/mu/vmin,-one*target_mass/(2*mu**2*vmin),-one*(delta*1e-6)**2/(2*target_mass*vmin)])
                    eft_nvel=[eft_amplitude_squared,eft_amplitude_squared1,eft_amplitude_squared1_e,eft_amplitude_squared1_em1]
                    for n_vel in n_vel_list:

                        diff_sigma_tilde=def_response_function_capture(coeff1,target,eft=lambda coeff1,q,element,n_isotope,coeff2: eft_modifier(q)*eft_nvel[n_vel](coeff1,q,element,n_isotope,coeff2,j_chi),coeff2=coeff2)

                        dc_dvi_de_cont=rho_loc/mchi*rho_target.reshape(-1,1)*(e2-e1)*p_vdep[n_vel]*diff_sigma_tilde(e1+t*(e2-e1))*lambda_i*hamiltonian.coeff_squared(__c1__,__c2__,**args_c1_c2)[tau,tau_prime]

                        if choose_zero:
                            dc_dvi_de_cont=np.choose(e2>e1,[0,dc_dvi_de_cont])
                        dc_dvi_de_target=dc_dvi_de_target+dc_dvi_de_cont 
        dc_dvi_de=dc_dvi_de+np.choose(w>vstar,[0,dc_dvi_de_target]) 
    return dc_dvi_de


def wimp_capture_accurate(celestial_body_obj,hamiltonian,vmin,delta_eta,mchi,j_chi=0.5,delta=0,rho_loc=0.3,response_functions=None,targets_list=None,coeff_squared_list=None,reset_response_functions=False,tau_range=range(2),tau_prime_range=range(2),n_vel_list=range(4),sum_over_streams=True,verbose=False,flatten=False,v_cut=0,choose_zero=True,**args): 
    '''
    Routine alternative to wimp_capture that calculates the WIMP capture rate without making use of the interpolation method.
    WARNING: significantly slower compared to wimp_capture 
    ---------------------------------------------------------
    Input: All its input arguments are in common with the wimp_capture routine
    ---------------------------------------------------------
    Output: The capture rate in s^{-1}.
    '''
    solar_radius_in_cm=6.9634e10 #to convert from solar radii to cm
    celestial_body_radius=celestial_body_obj.radius*solar_radius_in_cm

    if targets_list is None:
        targets_list=celestial_body_obj.targets

    func2=lambda t : d_wimp_capture_dvi_de(celestial_body_obj,hamiltonian,t,vmin,delta_eta,mchi,j_chi=j_chi,delta=delta,rho_loc=rho_loc,targets_list=targets_list,coeff_squared_list=coeff_squared_list,tau_range=tau_range, tau_prime_range=tau_prime_range,n_vel_list=n_vel_list,v_cut=v_cut,choose_zero=choose_zero,**args)

    sigma_int2,err=integrate.quad_vec(func2,0,1)

    wimp_capture_accurate.dc_dvi=sigma_int2

    cc_acc=celestial_body_radius**3*np.sum(integrate.simps(4*np.pi*(sigma_int2.T*celestial_body_obj.r_vec**2),celestial_body_obj.r_vec))

    return cc_acc


def e_intervals_capture(target,celestial_body,mchi,delta,vmin,v_cut=0):
    if vmin[0]==0:
        vmin=vmin[1:]

    mu=1./(1./mchi+1./target.mass)

    vstar=get_vstar(target.mass, mchi, delta)

    v_esc=celestial_body.v_esc
    w=np.sqrt(v_esc.reshape(-1,1)**2+vmin**2)
    w_larger_than_vstar=np.choose(w>vstar,[vstar,w])

    e1_temp=er_min(target.mass, mchi, delta, w_larger_than_vstar)

    ecap=er_capture(mchi,delta,vmin,v_cut=v_cut)

    e1=np.choose(e1_temp>ecap,[ecap,e1_temp])
    e2=er_max(target.mass, mchi, delta, w_larger_than_vstar)

    return e1,e2

def get_wimp_capture_response_functions_energy_sampling(target, celestial_body_list, mchi_list,delta_list,vmin, n_sampling=1000,**args): 

    '''
    calculates an optimized list of energy values in keV that can be used to increase the sampling of the integrated response functions sigma_tilde. See WimPyC paper(arXiv:2510.21185) for details of how the optimization is performed.
    ---------------------------------------------------------
    Input:

    - target              - object belonging to the isotope class
    
    - celestial_body_list - list or array of objects belonging to the celestial_body class

    - mchi_list           - list or array of values of the WIMP mass in GeV

    - delta_list          - list or array of value of the mass splitting delta in keV 
    
    - vmin                - list or array of WIMP asymptotic speeds (at large distance from the celestial body) in km/s

    - n_sampling          - (approximate) lenth of output array
    ---------------------------------------------------------
    Output:  array of energy values in keV

    '''



    er_sampling=np.array([])
    for celestial_body in celestial_body_list:
        for mchi in mchi_list:
            for delta in delta_list:
                if target.symbol in [t.symbol for t in celestial_body.targets]:
                    e1_temp,e2_temp=e_intervals_capture(target,celestial_body,mchi,delta,vmin)

                    e1=e1_temp[e1_temp<e2_temp]
                    e2=e2_temp[e1_temp<e2_temp]
                    er_sampling=np.append(er_sampling,np.unique(np.append(np.ndarray.flatten(e1),np.ndarray.flatten(e2))))
    if n_sampling==0:
        return np.sort(er_sampling)
    else:
        n_gap=max(int(len(er_sampling)/n_sampling),1)
        return np.sort(er_sampling)[::n_gap]

def merge_tuples(tuple1,tuple2):
    tuple_n_c1c2=()
    for n_c1c2 in range(len(tuple1)):
        tuple_n_vel=()
        for n_vel in range(len(tuple1[0])-1):
            tuple_tau=()
            for tau in range(len(tuple1[0][0])):
                tuple_tau_prime=()
                for tau_prime in range(len(tuple1[0][0][0])):
                    er1,sigma1=tuple1[n_c1c2][n_vel][tau][tau_prime]
                    er2,sigma2=tuple2[n_c1c2][n_vel][tau][tau_prime]

                    er12,nn=np.unique(np.append(er1,er2),return_index=True)
                    sigma12=np.append(sigma1,sigma2)[nn]

                    er12_sorted=np.sort(er12)
                    sigma12_sorted=sigma12[np.argsort(er12)]

                    tuple_tau_prime+=[er12_sorted,sigma12_sorted],
                tuple_tau+=tuple_tau_prime,
            tuple_n_vel+=tuple_tau,
        tuple_n_vel+=tuple1[n_c1c2][-1],
        tuple_n_c1c2+=tuple_n_vel,

    tuple12=tuple_n_c1c2
    return tuple12

def dump_response_functions_capture(isotope_obj,hamiltonian,j_chi,tuple,filename_suffix=''):
    path_response_functions=WIMP_CAPTURE_PATH+'/Response_functions'

    path_spin=path_response_functions+'/spin_'+print_spin(j_chi)
    path=path_spin

    for n_c1c2,(__c1__,__c2__) in enumerate(hamiltonian.coeff_squared_list):

        outputfile=isotope_obj.symbol+'_c_'+print_coupling(__c1__)+'_c_'+print_coupling(__c2__)+filename_suffix+'.npy'
        with open(path+'/'+outputfile, 'wb') as file:
            pickle.dump(tuple[n_c1c2], file)

def wimp_capture_geom(celestial_obj,mchi,vmin,delta_eta,rho_loc=0.3):
    '''
    calculates the geometrical capture rate
    ---------------------------------------------------------
    Input:
    -  celestial_obj      - object belonging to celestial_body class 
    -  mchi               - WIMP mass in GeV
    -  vmin               - an array containing a list of WIMP asymptotic speeds in km/s.
    -  delta_eta          - an array containing the contribution of each stream to the halo 
                             function eta(v) in (km/sec)^-1 (the routine streamed_halo_function can be used to calculate it) 
    - rho_loc             - local dark matter density in units of GeV/cm^3.
                            Default: 0.3

    ---------------------------------------------------------
    Output:  WIMP geometrical capture in sec^{-1}

    '''

    solar_radius_in_cm=6.9634e10
    R_star=celestial_obj.radius*solar_radius_in_cm*1e-5 #in km
    v_esc_star=celestial_obj.v_esc[-1] #in km/s
    focusing_effect=np.sum(vmin**2*delta_eta)+v_esc_star**2*np.sum(delta_eta)

    return np.pi*R_star**2*rho_loc/mchi*focusing_effect*1e15 #in s**-1

def wimp_capture_annihilation(celestial_body_obj, mchi, capture_rate, sigma_v=3e-26, t_age=4.603e9, effective_volume=None):
    '''
    calculates the WIMP annihilation rate following K. Griest and D. Seckel, Cosmic Asymmetry, Neutrinos and the Sun, Nucl. Phys. B 283 (1987) 681–705.
    ---------------------------------------------------------
    Input:
    - celestial_obj       - object belonging to celestial_body class 
    - capture_rate        - wimp capture rate in sec^{-1}
    - mchi                - WIMP mass in GeV
    - t_age               - age of the celestial body in years.
                             Default value: 4.603e9 (age of the Sun)
    - sigma_v             - WIMP annihilation cross section times velocity (lim v->0)
                             Default: 3e-26 cm^3 sec ^{-1}
    - effective_volume    - if None, the standard thermalization model of K. Griest and D. Seckel is used. For alternative scenarios, see WimPyC paper(arXiv:2510.21185). 
    ---------------------------------------------------------
    Output:  WIMP annihilation rate in sec^{-1}

    '''
    if  effective_volume is None:

        V1=6.5e28*(celestial_body_obj.T_c/1.4e7/(mchi/10))**(3/2) #cm^3

        V2=6.5e28*(celestial_body_obj.T_c/1.4e7/(2*mchi/10))**(3/2) #cm^3

        
        effective_volume=V1**2/V2

    C_A=sigma_v/effective_volume
    

    t_age=t_age*31536000 #conversion from years to seconds
    t_eq=1/np.sqrt(capture_rate*C_A) #second

    wimp_capture_annihilation.t_eq=t_eq
    
    return capture_rate/2*np.tanh(t_age/t_eq)**2





Sun=celestial_body('Sun')
Earth=celestial_body('Earth')
Jupiter=celestial_body('Jupiter')
White_Dwarf=celestial_body('White_Dwarf')
MS_Star=celestial_body('MS_Star')

celestial_body_list=[Sun,Earth,Jupiter,White_Dwarf,MS_Star]
celestial_body_obj_names_string='WD.'+' WD.'.join([cb.name for cb in celestial_body_list])
def list_celestial_bodies():
    print('Available celestial body objects in WimPyDD:')
    print('import WimPyDD as WD')
    print(celestial_body_obj_names_string)
    print('Use print() to get info on each celestial_body. For instance:')
    print('Type print(WD.Sun)')
    print(Sun)

##define all isotopes contained in celestial bodies
targets_list=np.array([])
for cb in celestial_body_list:
    targets_list=np.append(targets_list,cb.targets)

targets_list_symbols=np.unique([t.symbol for t in targets_list])
isotope_names_list=[''.join([t for t in symbol if t.isalpha()])+''.join([t for t in symbol if t.isdigit()]) for symbol in targets_list_symbols]

for name,symbol in zip(isotope_names_list,targets_list_symbols):
    exec(name+"=isotope('"+symbol+"')")


isotopes_list=np.array([])
for name,symbol in zip(isotope_names_list,targets_list_symbols):
    exec(name+"=isotope('"+symbol+"')")
    exec('isotopes_list=np.append(isotopes_list,'+name+')')
    
def list_isotopes():
    print('Available isotope objects in WimPyDD:')
    print('import WimPyDD as WD')
    print('WD.'+' WD.'.join(isotope_names_list))
    print('Use print() to get info on each isotope. For instance:')
    print('Type print(WD.Fe56)')
    print(Fe56)
#
