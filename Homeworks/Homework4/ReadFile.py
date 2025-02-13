import numpy as np
import astropy.units as u

def Read(filename):
    """
        Reads in particle data for model of Milky Way at a specific SnapNumber. 
        
        Inputs:
            filename (string) of where the data text file is stored, e.g. 'MW_000.txt'

        Outputs:
            time (astropy units Myr) corresponding to SnapNumber of model. 

            N_tot is the total number of particles across all types in the model. 

            data (numpy array) is a table describing the following for each particle in the modeL:
                - type: Dark Matter [labeled 1], Disk Stars [labeled 2], or Bulge Stars [labeled 3]
                - m: Mass of particle (astropy units of 10^10 Msun
                - x: x position (astropy units of kpc) of particle from center of mass of MW 
                - y: y position (astropy units of kpc) of particle from center of mass of MW 
                - z: z position (astropy units of kpc) of particle from center of mass of MW 
                - vx: x-component velocity (astropy units of km/s) of particle in Cartesian frame centered on the MW 
                - vy: y-compent velocity (astropy units of km/s) of particle in Cartesian frame centered on the MW 
                - vz: z-component velocity (astropy units of km/s) of particle in Cartesian frame centered on the MW 
    """

    file = open(filename, 'r') # Opens .txt file to access text
    line1 = file.readline()  # Reads in first line from file where the time information is stored
    label, value = line1.split() # Extracts the "Time" label and its corresponding value from line1
    time = float(value) * u.Myr # Converts time value from a string to a float in Myr

    line2 = file.readline() # Reads in second line from file where the number of particles information is stored
    label, value = line2.split() # Extracts the "Total" label and its corresponding value from line2
    N_tot = float(value) # Converts total number of particles value from a string to a float 

    file.close() # Closes file
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3) # Formats particle info table into a numpy array with the same column names, skipping the first 3 lines of .txt file
    return time, N_tot, data