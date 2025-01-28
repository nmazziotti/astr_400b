import numpy as np
import astropy.units as u

def ReadFile(filename):
    file = open(filename, 'r')
    line1 = file.readline()
    label, value = line1.split()
    time = float(value) * u.Myr

    line2 = file.readline()
    label, value = line2.split()
    N_tot = float(value) 

    file.close()
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    return time, N_tot, data