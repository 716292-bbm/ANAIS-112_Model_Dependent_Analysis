from Final.Librerias.GeneraRitmo.Libs.RAPIDD.rapidd.core import base_dir, reset_coefficients, isofromneuc, set_any_Ncoeff, read_halo, difrate_dER
import numpy as np

def calc_ANAIS_rates(dmmass, cp=1.0e-4, cn=1.0e-4, op=1, rhoDM=0.4, output_path="output", halo_path=base_dir+"/SHM.dat"):
    

    reset_coefficients()

#     rhoDM = card['rhoDM']  
    read_halo(halo_path)


    set_any_Ncoeff(cp, op, "p")
    set_any_Ncoeff(cn, op, "n")

    
    energies = np.linspace(0.01, 40, 100)

#     drde_Xe = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Xe", basis="All")
#     drde_Ar = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ar", basis="All")
#     drde_Ge = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ge", basis="All")

    drde_Na = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Na", basis="All") 
    drde_I = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="I", basis="All")


    np.savez(output_path+"/DDrates.npz", energies=energies, drde_Na=drde_Na, drde_I=drde_I)

    np.savetxt(output_path+"/DDrates.txt", np.column_stack((energies,drde_Na, drde_I)),
           header="Er drde_Na drde_I")

    return


def plot_ANAIS_rates(dmmass, cp, cn, op, rhoDM=0.4, halo_path=base_dir+"/SHM.dat"):
    

    import matplotlib.pyplot as plt
    
    reset_coefficients()

#     rhoDM = card['rhoDM']  
    read_halo(halo_path)


    set_any_Ncoeff(cp, op, "p")
    set_any_Ncoeff(cn, op, "n")

    
    energies = np.linspace(0.01, 40, 100)

#     drde_Xe = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Xe", basis="All")
#     drde_Ar = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ar", basis="All")
#     drde_Ge = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Ge", basis="All")

    drde_Na = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Na", basis="All") 
    drde_I = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="I", basis="All")

    plt.loglog(energies, drde_Na, label='Na')
    plt.loglog(energies, drde_I, label='I')

    plt.ylabel(r'$dR/dE_{R}\,\,\left[{\rm kg}\,{\rm day}\,{\rm keV}\right]^{-1}$') 
    plt.xlabel(r'$E_{R}\,\left[{\rm keV}\right]$')
    plt.legend()
    plt.show()
    return

def plot_ANAIS_multiop(dmmass, cp, cn, rhoDM=0.4, 
                       halo_path=base_dir + "/SHM.dat",
                       ncols=3, figsize=(15,10)):
    """
    Plot ANAIS differential rates for multiple NR operators in subplots.
    
    Parameters
    ----------
    dmmass : float
        Dark matter mass in GeV.
    cp, cn : float
        Proton and neutron coefficients.
    op_list : list of int
        List of operators to loop over, e.g. [1, 5, 8, 11].
    rhoDM : float
        Local DM density.
    halo_path : str
        Path to halo file.
    ncols : int
        Number of subplot columns.
    figsize : tuple
        Figure size.
    """
    
    import matplotlib.pyplot as plt

    
    # Prepare figure
    op_list = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    nops = len(op_list)
    nrows = int(np.ceil(nops / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    
    energies = np.linspace(0.01, 40, 200)

    # Load halo once
    reset_coefficients()
    read_halo(halo_path)

    # Loop over operators
    for i, op in enumerate(op_list):
        ax = axes[i // ncols][i % ncols]

        # Reset coefficients each subplot
        reset_coefficients()
        set_any_Ncoeff(cp, op, "p")
        set_any_Ncoeff(cn, op, "n")

        # Compute rates
        drde_Na = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="Na", basis="All")
        drde_I  = np.vectorize(difrate_dER)(rhoDM, dmmass, energies, target="I",  basis="All")

        # Plot
        ax.loglog(energies, drde_Na, label="Na")
        ax.loglog(energies, drde_I,  label="I")

        ax.set_title(f"Operator O{op}")
        ax.set_xlabel(r"$E_R$ [keV]")
        ax.set_ylabel(r"$dR/dE_R$ [kg day keV]$^{-1}$")
        ax.legend(fontsize=9)

    # Remove empty panels
    for j in range(nops, nrows * ncols):
        fig.delaxes(axes[j // ncols][j % ncols])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_ANAIS_multiop(dmmass=50.0, cp=1e-4, cn=1e-4)

    