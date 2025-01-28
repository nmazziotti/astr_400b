import numpy as np
import astropy.units as u

from ReadFile import ReadFile

def ParticleInfo(filename, ParticleType, ParticleNum):
    
    time, N_tot, data = ReadFile(filename)

    if ParticleType == 'DM':
        index = np.where(data['type'] == 1.0)
    elif ParticleType == 'Disk':
        index = np.where(data['type'] == 2.0)
    elif ParticleType == 'Bulge':
        index = np.where(data['type'] == 3.0)

    if ParticleNum >= len(data[index]) or ParticleNum < 0:
        raise Exception(f'Invalid particle number for this particle type.')
    else:
        Particle = data[index][ParticleNum]
        mass = Particle['m'] * 1e10 * u.M_sun
        position = Particle[['x', 'y', 'z']]
        velocity = Particle[['vx', 'vy', 'vz']]
        distance = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2) * u.kpc
        speed = np.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2) * u.km/u.s

        return np.round(distance,3), np.round(speed, 3), mass