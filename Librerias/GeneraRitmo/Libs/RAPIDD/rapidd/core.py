import ctypes

import numpy as np
import os
import platform


# Determine the platform and set the library name accordingly
if platform.system() == 'Linux':
    lib_name = 'libRAPIDD.so'
elif platform.system() == 'Darwin':  # Darwin is macOS
    lib_name = 'libRAPIDD.dylib'
elif platform.system() == 'Windows':
    lib_name = 'RAPIDD.dll'
else:
    raise OSError("Unsupported operating system")

#### Call the shared library #######

base_dir = os.path.dirname(os.path.realpath(__file__))

# Construct the absolute path to the shared library
lib_path = os.path.join(base_dir, '..', 'lib', 'build', lib_name)


_crapidd = ctypes.CDLL(lib_path)


vev = 246.2
mneutron = 0.939565
mproton = 0.938272



########## Halo related stuff ################

_define_and_write_halo_path = _crapidd.define_and_write_halo_path

_define_and_write_halo_path.argtypes = [ctypes.c_char_p, ctypes.c_char_p,ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double ,ctypes.c_double, ctypes.c_double ,ctypes.c_int ]

_define_and_write_halo_path.restype = ctypes.c_void_p

### Vsolar = v0 + vpec which we put that to ve ignoring earths orbit for now 
#  https://arxiv.org/pdf/2105.00599 use table 

vearth_ref = np.sqrt(11.1**2 + (238.0+12.2)**2 + 7.3**2) 

def calc_new_halo(path, profile="SHM", vesc=544, v0=238, beta=0, vt=90., vc=245, ve=vearth_ref, k=0.0, i=2 ):
    _define_and_write_halo_path(path.encode(), profile.encode(), vesc, v0, beta, vt, vc, ve, k, i)
    return

_read_halo = _crapidd.read_halo

_read_halo.argtypes = [ctypes.c_char_p]

_read_halo.restype = ctypes.c_void_p 

halo_path = os.path.join(base_dir, '../lib/halo_table.dat')


def read_halo(path = halo_path):
    _read_halo(path.encode())
    return


###### Setting coefficients ######




_set_any_coeffs = _crapidd.set_any_coeffs

_set_any_coeffs.argtypes = [ctypes.c_double, ctypes.c_int]

def set_any_coeffs(C, i):
    return _set_any_coeffs(C, i)


_set_any_Ncoeff = _crapidd.set_any_Ncoeff
_set_any_Ncoeff.argtypes= [ ctypes.c_double, ctypes.c_int, ctypes.c_char_p]
_set_any_Ncoeff.restype = ctypes.c_void_p


_Cp = _crapidd.Cp
_Cp.argtypes = [ctypes.c_int]
_Cp.restype = ctypes.c_double


_Cn = _crapidd.Cn
_Cn.argtypes = [ctypes.c_int]
_Cn.restype = ctypes.c_double

def set_any_Ncoeff(coeff, op, nuc):
    '''coeff, op, nuc'''
    _set_any_Ncoeff(coeff, op, nuc.encode())
    return 

def reset_coefficients() :
    for i in range(16):
        set_any_Ncoeff(0, i, "p")
        set_any_Ncoeff(0, i, "n")
    
    return

def Cp_val(op):
    return _Cp(op)

def Cn_val(op):
    return _Cn(op)


def isofromneuc(cp, cn):
    ''' Simply takes coeffs from p n basis to 0 1 basis '''
    return (cp+cn)/2, (cp-cn)/2

####### differential rate #########
_difrate_dER_python = _crapidd.difrate_dER_python
_difrate_dER_python.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
_difrate_dER_python.restype = ctypes.c_double

@np.vectorize
def difrate_dER(rhoDM, mDM, ER, model="None", target="Xe", basis="ISO"):
    return _difrate_dER_python(rhoDM, mDM, np.log10(ER), model.encode(), target.encode(), basis.encode())

