from numpy import pi, sqrt, exp, piecewise, linspace
from scipy.special import erf



def eta_shm(vmin, v0, vesc, vE):

    '''analytic SHM eta arXiv:1509.01598'''

    K = v0**3 * pi * (sqrt(pi) * erf(vesc/v0) - 2 * (vesc/v0)* exp(-(vesc/v0)**2))

    def eta1(vmin): return ((v0**2 *pi )/(2*vE*K)) * ( -4 * vE * exp(-(vesc/v0)**2) + sqrt(pi) * v0 *( erf((vmin + vE)/v0) - erf((vmin-vE)/v0) ) )
    def eta2(vmin): return ((v0**2 *pi )/(2*vE*K)) * (  -2 * (vesc - vmin + vE) * exp(-(vesc/v0)**2) + sqrt(pi) * v0*( erf(vE/v0) - erf((vmin-vE)/v0) )  )

    conds = [vmin < vesc - vE , (vesc - vE < vmin) & (vmin < vesc + vE)]

    return piecewise(vmin, conds, [eta1, eta2])


def gen_shm_table(path, card):
    '''Generates the SHM eta function through analytic expression'''

    # --- Read this from maddm_card ---
    v0 = card['vmp'] # km/s
    vesc = card['vescape'] # km/s 

    # --- Read this from updated maddm_card --- 

    #vE_HC = [29.2, -0.1, 5.9] # Earth velocity relative to the Sun (heliocentric) at March 9 in km/s
    #vS = [11.1, 12.2, 7.3] # Solar peculiar velocity in km/s

    vE_HC = card['vEarth_mod']*[0.97986577, - 0.0033557 ,  0.19798658] # Earth velocity relative to the Sun (heliocentric) at March 9 in km/s 
    
    vS = [card['vSun_r'], card['vSun_phi'], card['cSun_theta']] # Solar peculiar velocity in km/s

    v0_vec = [0, v0, 0]

    # Calculate the module of Earth velocity relative to the galactic center
    vE_GC = sqrt(
        (vE_HC[0] + vS[0] + v0_vec[0])**2 +
        (vE_HC[1] + vS[1] + v0_vec[1])**2 +
        (vE_HC[2] + vS[2] + v0_vec[2])**2
    )

    vend = vesc + vE_GC
    vminspace = linspace(0.0, vend, 100)
    eta_values = eta_shm(vminspace, vE_GC, v0, vesc)


    # Open the file for writing
    with open(path, "w") as f:
        # First line: number of entries and a "1"
        f.write(f"{len(vminspace)} 0\n")
        
        # Write vmin and eta values line by line
        for v, eta in zip(vminspace, eta_values):
            f.write(f"{v:.5E} {eta:.5E}\n")