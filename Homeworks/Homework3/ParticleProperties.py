import numpy as np
import astropy.units as u

from ReadFile import ReadFile

def ParticleInfo(filename, ParticleType, ParticleNum):
    """
        Computes mass, distance, and velocity of a specific particle within a given type (DM, Disk, Bulge).

    Inputs:
        filename (string) of where the data text file is stored, e.g. 'MW_000.txt'

        ParticleType (string): component of the MW that a particle is associated with, either 
            DM, Disk, or Bulge.

        ParticleNum (int): Index of particle within the given particle type, NOT index of the particle in the table. 

    Outputs:
        distance (astropy units kpc) of particle from center of mass of MW

        velocity (astropy units km/s): magnitude of the particle's 3D velocity in Cartesian frame centered on the MW.

        mass (astropy units 1e10 Msun) of selected particle 
    """
    
    time, N_tot, data = ReadFile(filename) # Extract model time, total number of particles, and data table of particles from .txt file

    # Filter data table to given part type (DM represented by 1.0, Disk is 2.0, and Bulge is 3.0)
    if ParticleType == 'DM':
        index = np.where(data['type'] == 1.0) 
    elif ParticleType == 'Disk':
        index = np.where(data['type'] == 2.0)
    elif ParticleType == 'Bulge':
        index = np.where(data['type'] == 3.0)

    # Raise exception if the given particle number does not exist 
    if ParticleNum >= len(data[index]) or ParticleNum < 0:
        raise Exception(f'Invalid particle number for this particle type.')
    else:
    # If particle number is valid, proceed as follows

        Particle = data[index][ParticleNum] # Extract row for this particle from the data table 
        mass = Particle['m'] * 1e10 * u.M_sun # Compute mass 
        pos_vec = Particle[['x', 'y', 'z']] # Store position coordinates as 3D vector
        v_vec = Particle[['vx', 'vy', 'vz']] # Store velocity components as 3D vector
        distance = np.sqrt(pos_vec[0]**2 + pos_vec[1]**2 + pos_vec[2]**2) * u.kpc # Compute distance from position vector
        velocity = np.sqrt(v_vec[0]**2 + v_vec[1]**2 + v_vec[2]**2) * u.km/u.s # Compute magnitude from velocity vector

        return np.round(distance,3), np.round(velocity, 3), mass