import uproot
import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy import special
from scipy.stats import norm


def carga_datos_hist():
    file = uproot.open("/home/borja/Documents/TFG/Datos/BEhistos_year123456.root")
    file.keys()
    longitud=290
    es=np.arange(1,30,0.1)
    datos=np.zeros((19,longitud))
    for i in range (9):
        datos[2*i+1]=file["hbea_123456y_D"+str(i)].values()
        datos[2*i+2]=file["hbea_123456y_D"+str(i)].errors()
    datos[0]=es
    return datos
    #np.savetxt("Exposicion_exp.csv",datos)   

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

#################
# Funcion inversa de la CDF de la distribucion gaussiana. (Devuelve el valor del limite superior de la integral gaussiana desde 
# menos infinito que da como resultado un area p (en tanto por uno),devuelve el valor en multiplos de sigma)
#################
# p: area de la gausiana para la cual calculamos el limite superior
def f_significancia(p):
    return norm.ppf(p)

#################
# Calcula el ritmo experimental -> Resultado en (c/kg/d)
#################
# E_i: Limite inferior de la energia (keV_ee) (Minimo 1 keV_ee)
# E_f: Limite superior de la energia (keV_ee) (Maximo 30 keV_ee)
# datos: Array con los datos de los histogramas en formato
#   Eje_y_H0 Error_y_H0 Eje_y_H1 Error_y_H1 ... Eje_y_H8 Error_y_H8 Eje_x

def integral_ritmo_exp(E_i,E_f,datos):

    if(E_i<1): 
        print('El valor de E_i debe ser mayor que 1 keV_ee')
        return "El valor de E_i debe ser mayor que 1 keV_ee"
    if(E_f>30): 
        print('El valor de E_f debe ser menor que 30 keV_ee')
        return "El valor de E_f debe ser menor que 30 keV_ee"
    

    #   Hacemos la integral en el intervalo de energias definido por E_i y E_f

    ritmos_exp=np.zeros(9)

    #   Iteramos en los 9 detectores
    for i in range (9):
        ritmos_exp[i]=0

        #   Sumamos en cada detector

        for j in range(int(10*(E_i-1)),int(10*(E_f-1))):
                ritmos_exp[i]+=datos[2*i+1][j]

        #   Multiplicamos por ancho de bin

        ritmos_exp[i]*=0.1

    return ritmos_exp
        
#################
# Calcula el numero de cuentas -> Resultado en (c)
#################
# Ritmo_exp: Array de 9 componentes con el ritmo en cada detector (c/kg/d)
# texp: Array de 9 componentes con el tiempo de exposicion en cada detector (d)

def numero_cuentas_exp(ritmo_exp,texp,array_det=np.ones(9,dtype=int)):
    retval=0
    for i in range (9):
        if(array_det[i]==1):
            retval+=ritmo_exp[i]*texp[i]*12.5
    return retval

#################
# Calcula el numero de cuentas -> Resultado en (c)
#################
# E_i: Limite inferior de la energia (keV_ee) (Minimo 1 keV_ee)
# E_f: Limite superior de la energia (keV_ee) (Maximo 30 keV_ee)
# texp: Array de 9 componentes con el tiempo de exposicion en cada detector (d)
# cl: Confidence Level

def numero_cuentas_exp_CL(E_i,E_f,cl=0.9,array_det=np.ones(9,dtype=int)):

    datos=carga_datos_hist()

    exposiciones=carga_exposiciones()

    #print(exposiciones)

    ritmos_exp=integral_ritmo_exp(E_i,E_f,datos)

    #print(ritmos_exp)

    numero_cuentas=numero_cuentas_exp(ritmos_exp,exposiciones,array_det)

    #print(numero_cuentas)

    k=f_significancia(cl)

    #print(k)

    retval=numero_cuentas+k*m.sqrt(numero_cuentas)

    return retval






