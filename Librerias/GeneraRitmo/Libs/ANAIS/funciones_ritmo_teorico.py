
import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy import special
from scipy.stats import norm

#ctes

mn=0.9315               # Masa nucleon (GeV/c^2)
NA=6.022e23             # Numero de Avogadro
c=3e8                   # Velocidad luz (m/s)
hbar = 6.582e-16        # eV*s
sqrtpi=m.sqrt(m.pi)     # Raíz de pi

#Halo DM
rho=0.3                 # Densidad local de materia oscura Gev/c2/cm3
vesc=544                # Velocidad de escape de la galaxia (km/s)

#Velocidad del sol
u1=np.array([ 0.9941, 0.1088, 0.0042])
u2=np.array([-0.0504, 0.4946, -0.8677])
v0=np.array([0,238,0])      
v0M=m.sqrt(v0.dot(v0))                      # Velocidad en reposo local standar (km/s)
v0n=np.sqrt(v0.dot(v0))                     # Norma de la Velocidad en reposo local standar (km/s)
omega=0.0172                                # Frecuencia angular en d-1 (2pi/365)
vpec=np.array([11.1,12.2,7.3])              # Velocidad peculiar (km/s)
vsun=v0+vpec                                # Velocidad del Sol (km/s)
vorb=29.8                                   # Velocidad orbital de la Tierra (km/s)
timeMax=72.9                                # Dia de la máxima velocidad de la tierra: 2 de junio, empezando desde el 22 de marzo.
timeMin=255.5                               # Dia de la mínima velocidad de la tierra: empezando desde el 22 de marzo.
timeAvg=(timeMax+timeMin)/2.                # Dia promedio (El que usaremos para calcular el ritmo)
ene=np.arange(0,100,1)                      # Array de energías de 0 a 100 MeV en intervalos de 1 MeV
time=np.arange(0,365,1)                     # Array de dias, 365 dias en intervalos de 1 dia
unitsvearth = 0 # km/s -> 0 , m/s -> 1


#################
# Funcion de Bessel de primera especie
#################
def j1(x):
    return(m.sin(x)/(x*x)-m.cos(x)/x)

#################
# Velocidad de la tierra, en km/s
#################
def vearth(t): 
    v = vsun + vorb*(m.cos(omega*t)*u1+m.sin(omega*t)*u2)
    return m.sqrt(v.dot(v))
    
#################
# Velocidad minima de la DM para un umbral E (en km/s)
# E en keV
# A: Numero másico
# mW: Masa del Wimp en GeV/c^2
#################
def vmin(E,A,mW): 
    mN=A*mn
    mu_N=mW*mN/(mW+mN)
    return m.sqrt(E*mN*1e-6/2./mu_N/mu_N)*c*1e-3 # km/s

#################
# mean inverse speed function (Función de velocidad inversa media) (en s/km)
# E en keV
# t: Tiempo en dias desde el 22 de marzo
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
#################
def eta(E,t,A,mW):
  x=vmin(E,A,mW)/v0M
  y=vearth(t)/v0M
  z=vesc/v0M
  N=m.erf(z)-2*z/sqrtpi*m.exp(-z*z)
  retval=1./2/y/v0M/N
  if x > z+y:
    return 0
  if x > z-y:
    retval*=m.erf(z)-m.erf(x-y)-2/sqrtpi*(z+y-x)*m.exp(-z*z)
  else:
    retval*=m.erf(x+y)-m.erf(x-y)-4/sqrtpi*y*m.exp(-z*z)
  return retval
        
#################
# Factor de Forma
# E en keV
# A: Numero másico
#################
def FF(E,A):
  if E==0:
    return 1
  R=1.2*A**0.3333 # fm
  s=1 # fm
  mN=A*mn
  R1=m.sqrt(R*R-5*s*s) # fm
  R1*=1e-6/hbar/c

  q=m.sqrt(2*1e-6*mN*E) # GeV
  x=q*R1
  #print("mN = " + str(mN) + " E= "+str(E)+ " q="+str(q)+" R1="+str(R1)+" x= " + str(x) + "hbarc= " +str(hbar*c*1e6))
  retval=m.sin(x)/x/x-m.cos(x)/x
  retval*=3/x
  retval*=retval
  aux=q*s*1e-6/hbar/c
  retval*=m.exp(-aux*aux)
  return retval

#################
# Ritmo diferencial  [cts/KeV/d/kg]
# E en keV
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
#################
def rate(E,t,A,mW,sigmaSI):
    Mdet=1000*NA*mn                                 # masa de 1 kg de detector [GeV/c^2]
    mu_n=mW*mn/(mW+mn)                              # masa reducida del sistema nucleon-WIMP [Gev/c^2]
    retval=Mdet*rho/2./mW/mu_n/mu_n*A*A*sigmaSI     # Producto de los primeros terminos [cts*c^2/Gev/cm/kg]
    retval*=c*c                                     # Conversion a 100*[cts*m/s^2/Gev/kg]
    retval*=FF(E,A)                                 # Multiplicamos por factor de forma atomico (Adimensional) 100*[cts*m/s^2/Gev/kg]
    retval*=eta(E,t,A,mW)                           # Multiplicamos por integral 0.1*[cts/GeV/s/kg]
    retval*=8.64e-3                                 # Convertimos a [cts/KeV/d/kg]
    
    return retval

#########################
#########################
# Ritmo total, integrado entre Ei y Ef (en c/kg/d)
# Ei, Ef: Energias inicial y final en keV
# t: tiempo en dias desde el 22 de marzo
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
def totalRate(Ei,Ef,t,A,mW,sigmaSI):
  ene=np.arange(Ei,Ef,0.1)
  rates=np.array([rate(e,t,A,mW,sigmaSI) for e in ene])
  return rates.sum()*0.1

#########################
#########################
# Ritmo en funcon del tiempo
# Ei, Ef: Energias inicial y final en keV
# t: tiempo en dias desde el 22 de marzo
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
def ratevsTime(Ei,Ef,A,mW,sigmaSI):
    dias = np.arange(0,365)
    ritmos_dias=([totalRate(Ei,Ef,t,A,mW,sigmaSI) for t in dias])
    return ritmos_dias



# FUNCIONES ESPECIFICAS PARA NAI

#########################
#########################
# Ritmo diferencial
# E en keV
# t: tiempo en dias desde el 22 de marzo
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
def RateNaI(E,t,mW,sigmaSI):
    ratesNa=rate(E,t,23,mW,sigmaSI)
    ratesI=rate(E,t,127,mW,sigmaSI)
    return (23.*ratesNa+127.*ratesI)/(23.+127.)

#########################
#########################
# Ritmo total, integrado entre Ei y Ef (en c/kg/d)
# Ei, Ef: Energias inicial y final en keV
# t: tiempo en dias desde el 22 de marzo
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
def totalRateNaI(Ei,Ef,t,mW,sigmaSI):
  ene=np.arange(Ei,Ef,0.1)
  rates=np.array([RateNaI(e,t,mW,sigmaSI) for e in ene])
  return rates.sum()*0.1

#########################
#########################
# Ritmo en funcon del tiempo
# Ei, Ef: Energias inicial y final en keV
# t: tiempo en dias desde el 22 de marzo
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
def ratevsTimeNaI(Ei,Ef,mW,sigmaSI):
    dias = np.arange(0,365)
    ritmos_dias=([totalRateNaI(Ei,Ef,t,mW,sigmaSI) for t in dias])
    return ritmos_dias
       




# FUNCIONES TENIENDO EN CUENTA EL FACTOR QUENCHING

#FACTOR QUENCHING
#   Devuelve el factor quenching del NaI evaluado en electron-equivalent energy ee
#   Parametros
# ee : float o array-like de Electron-equivalent energy.
# Devuelve: float or np.ndarray
# Quenching factor QF(ee), interpolado linealmente en (x=ER*QF, y=QF)
# Construido para ER in [1, 100] with N=200 points.
# From modified Lindhard Fit of Na QF data (Tamara's thesis)
# k = 0.072, alpha = 0.007
# epsilon = alpha * ER
# g(ER)=3*pow(epsilon, 0.15) + 0.7*pow(epsilon,0.6) + epsilon
# QF(ER)=k*g(ER)/(1+k*g(ER))
# Eee=ER*QF(ER)

def getQFNa(ee, k=0.072, alpha=0.007, N=200, ER_min=1.0, ER_max=100.0):
    # calculate QF(ER) for interval ER_min, ER_max
    ER = np.linspace(ER_min, ER_max, N)
    epsilon = alpha * ER
    g = 3.0 * epsilon**0.15 + 0.7 * epsilon**0.6 + epsilon
    qf = (k * g) / (1.0 + k * g)

    # Graph: x = ER * qf (== ee), y = qf
    x = ER * qf
    y = qf

    # Interpolate y at the provided ee values.
    # For values outside the tabulated range, clamp to the nearest endpoint
    ee_arr = np.atleast_1d(ee).astype(float)
    y_interp = np.interp(ee_arr, x, y, left=y[0], right=y[-1])

    return y_interp[0] if np.isscalar(ee) else y_interp

def getQFI(ee, p0=0.03, p1=0.0006, limit=80.0, N=200, ER_min=1.0, ER_max=100.0):
    # calculate QF(ER) for interval ER_min, ER_max
    ER = np.linspace(ER_min, ER_max, N)

    # Piecewise-linear QF(ER)
    qf = p0 + p1 * np.minimum(ER, limit)

    # Graph x = ER * QF(ER), y = QF(ER)
    x = ER * qf
    y = qf

    # Interpolate y at ee; clamp outside domain to the nearest endpoint
    ee_arr = np.atleast_1d(ee).astype(float)
    y_interp = np.interp(ee_arr, x, y, left=y[0], right=y[-1])

    return y_interp[0] if np.isscalar(ee) else y_interp

#################
# Ritmo diferencial (Teniendo en cuenta el Quenching, Energía en ee) (en c/kevee/kg/d) 
# E en keV
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
#################
def rate_ee(Eee,t,A,mW,sigmaSI,Q=1):
    E=Eee/Q
    retval=rate(E,t,A,mW,sigmaSI)
    retval/=Q
    return retval

#########################
#########################
# Ritmo diferencial en detector NaI (en c/kevee/kg/d) 
# E: en keVee
# t: tiempo en dias desde el 22 de marzo
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
# QNa es el factor de Quenching a energía Eee(keVee) para Na
# QI es el factor de Quenching a energía Eee(keVee) para I
def rateNaI_ee(Eee,t,mW,sigmaSI, QNa=1, QI=1):
  # convert Eee to keV_NR
  rateNa = rate_ee(Eee,t,23,mW,sigmaSI,QNa)
  rateI = rate_ee(Eee,t,127,mW,sigmaSI,QI)
 
  return (rateNa*23.+rateI*127.)/(23.+127.)

#########################
#########################
# Ritmo total, integrado entre Eiee y Efee (en c/kg/d)
# Eiee, Efee: Energías inicial y final en keVee
# t: tiempo en dias desde el 22 de marzo
# A: numero masico
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2

def totalRate_NaI_ee(Eiee,Efee,t,mW,sigmaSI):
  energy_ee = np.arange(Eiee,Efee,0.1)
  qfNa = getQFNa(energy_ee)
  qfI = getQFI(energy_ee)
  rates = np.array([
        rateNaI_ee(e, t, mW, sigmaSI, qna, qi)
        for e, qna, qi in zip(energy_ee, qfNa, qfI)
    ])

  return rates.sum()*0.1

#########################
#########################
# Ritmo total, integrado entre Eiee y Efee (en c/kg/d)
# Eiee, Efee: Energías inicial y final en keVee
# t: tiempo en dias desde el 22 de marzo
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2

def totalRate_NaI_ee_DAMA(Eiee,Efee,t,mW,sigmaSI):
  energy_ee = np.arange(Eiee,Efee,0.1)
  qfNa = 0.3
  qfI = 0.09
  rates = np.array([
        rateNaI_ee(e, t, mW, sigmaSI, qfNa, qfI)
        for e in energy_ee
    ])

  return rates.sum()*0.1

def carga_exposiciones():
    exposiciones=np.zeros(9) #dias
        
    exposiciones[0]=2031.38
    exposiciones[1]=2033.20
    exposiciones[2]=2029.52
    exposiciones[3]=2022.55
    exposiciones[4]=2033.01
    exposiciones[5]=2030.18
    exposiciones[6]=2032.27 
    exposiciones[7]=2031.02
    exposiciones[8]=2020.29

    return exposiciones
def calcula_t_exposicion(det): #  Calculamos el tiempo total de exposicion (dias)

    exposiciones=carga_exposiciones()
    
    return exposiciones[det]

def calcula_m_exposicion(det): #  Calculamos la masa total de exposicion (kg)        
    return 12.5

#########################
#########################
# Numero de cuentas total, integrado entre Eiee y Efee (en cts)
# Eiee, Efee: Energías inicial y final en keVee
# t: tiempo en dias desde el 22 de marzo
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
# texp: Tiempo de exposicion en dias
# mexp: Masa del detector de NaI en kg

def numero_cuentas_teo(Eiee,Efee,t,mW,sigmaSI,array_det=np.ones(9,dtype=int)):
  Matriz=np.eye(9)
  exp=0
  for i in range(9):
          if(array_det[i]==1):
            texp=calcula_t_exposicion(i)
            mexp=calcula_m_exposicion(i)
            exp+=texp*mexp
  retval=totalRate_NaI_ee(Eiee,Efee,t,mW,sigmaSI)
  return retval*exp

#########################
#########################
# Numero de cuentas total, integrado entre Eiee y Efee (en cts)
# Eiee, Efee: Energías inicial y final en keVee
# t: tiempo en dias desde el 22 de marzo
# mW: Masa del Wimp en GeV/c^2
# sigmaSI: Seccion eficaz spin independent en cm^2
# texp: Tiempo de exposicion en dias
# mexp: Masa del detector de NaI en kg

def numero_cuentas_teo_DAMA(Eiee,Efee,t,mW,sigmaSI,array_det=np.ones(9,dtype=int)):
  Matriz=np.eye(9)
  exp=0
  for i in range(9):
          if(array_det[i]==1):
            texp=calcula_t_exposicion(i)
            mexp=calcula_m_exposicion(i)
            exp+=texp*mexp
  retval=totalRate_NaI_ee_DAMA(Eiee,Efee,t,mW,sigmaSI)
  return retval*exp