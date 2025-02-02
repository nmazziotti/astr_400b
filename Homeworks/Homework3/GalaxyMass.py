import numpy as np
from ReadFile import Read

def ComponentMass(filename, ParticleType):
    """
    Sums over all particles of the given galaxy component to compute the total mass.

    Inputs:
        filename (string) corresponding to the text file of the galaxy model

        ParticleType (string) indicating the desired mass component (DM, Disk, or Bulge)

    Outputs:
        M_tot (astropy units) is the total mass of the galaxy component in 1e12 Msun. 
    """

    time, N_tot, data = Read(filename)

    # Filter data table to given part type (DM represented by 1.0, Disk is 2.0, and Bulge is 3.0)
    if ParticleType == 'DM':
        index = np.where(data['type'] == 1.0) 
    elif ParticleType == 'Disk':
        index = np.where(data['type'] == 2.0)
    elif ParticleType == 'Bulge':
        index = np.where(data['type'] == 3.0)

    M_tot = np.sum(data[index]['m']) / 100

    return np.round(M_tot, 3)
