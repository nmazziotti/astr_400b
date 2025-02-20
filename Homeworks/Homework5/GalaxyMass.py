import numpy as np
from ReadFile import Read

def ComponentMass(filename, ParticleType):
    """
    Sums over all particles of the given galaxy component to compute the total mass.

    Inputs:
        filename (string) corresponding to the text file of the galaxy model

        ParticleType (int) indicating the desired mass component (DM = 1, Disk = 2, or Bulge = 3)

    Outputs:
        M_tot (astropy units) is the total mass of the galaxy component in 1e12 Msun. 
    """

    time, N_tot, data = Read(filename) # Retrieve time of SnapNumber, total number of particles in model, and particle data table

    # Filter data table to given part type (DM represented by 1.0, Disk is 2.0, and Bulge is 3.0)
    index = np.where(data['type'] == ParticleType) 
    
    # Sum over all particle masses in the array to get total mass for the component
    M_tot = np.sum(data[index]['m']) / 100 # Data is already stored in units of 1e10 Msun, so divide by 100 to get 1e12 Msun

    return np.round(M_tot, 3) # Round to 3 decimal places
