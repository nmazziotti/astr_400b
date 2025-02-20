import numpy as np 
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.constants import G

from CenterOfMass import CenterOfMass
from ReadFile import Read

# Convert G to convenient units of (kpc * km^2)/(s^2 & Msun)
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

class MassProfile:
    
    def __init__(self, galaxy, snap):
        """
        Class to calculate the radial mass and circular velocity profiles of a galaxy 
        according to a specific mass component or the total mass. 
        
            PARAMETERS
            ----------
            galaxy (string): Notation of galaxy (e.g. MW, M31, or M33)
            snap (int): Snapshot number in simulation (0 = present day) 
        """

        # Construct filename for simulation data using galaxy and snap
        # add a string of the filenumber to the value “000”
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename="%s_"%(galaxy) + ilbl + '.txt'
        
        # Store name of galaxy as gname 
        self.gname = galaxy

        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)                                                                                             
        # store the mass and positions of all particles 
        # mass in file is stored in units of 1e10 Msun
        self.m = self.data['m'] * 1e10 # leave out units of Msun for computing purposes 
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc
         
    
    def MassEnclosed(self, Radius, ptype):
        """
        Method that computes the mass enclosed within a set a radial distances from a galaxy's COM position.
        Only considers the mass of a given particle type. 
        
        INPUTS:
            Radius (numpy array of astropy quantities): Array of distances from the COM position in astropy units of kpc
            ptype (int): Integer corresponding to a particle type (1 = Halo, 2 = Disk, 3 = Bulge)
            
        OUTPUTS:
            masses (numpy array of astropy quantities): Array of the mass enclosed for each corresponding distance in 
            Radius array. Astropy units of Msun. 
        """
        
        COM = CenterOfMass(self.filename, 2) # construct CenterOfMass object using disk stars 
        COM_P = COM.COM_P() # compute COM position using COM_P method
        
        # select only particles of the given ptype and their masses 
        index = np.where(self.data['type'] == ptype)
        m_new = self.m[index] 
        
        # recalculate new particle positions in COM frame 
        x_new = self.x[index] - COM_P[0]
        y_new = self.y[index] - COM_P[1]
        z_new = self.z[index] - COM_P[2]
        
        # Calculate distance of each particle from the COM position in kpc. 
        COM_dists = np.sqrt(x_new**2 + y_new**2 + z_new**2)
        
        # Construct empty numpy array to store mass enclosed at each distance in Radius
        masses = np.zeros(len(Radius))
        
        # Loop over each distance in Radius array 
        for i in range(len(Radius)):
            enclosed = np.where(COM_dists < Radius[i]) # select only particles enclosed by this radius 
            mass_enclosed = np.sum(m_new[enclosed]) # Sum particle masses to calculate mass enclosed 
            masses[i] = mass_enclosed # Store mass enclosed in masses array 
        
        # Return masses array in astropy units of Msun 
        return np.array(masses) * u.Msun
    
    def TotalMassEnclosed(self, Radius):
        """
        Method that computes the total mass enclosed within a set a radial distances from a galaxy's COM position.
        
        INPUTS:
            Radius (numpy array of astropy quantities): Array of distances from the COM position in astropy units of kpc.
            
        OUTPUTS:
            TotalMasses (numpy array of astropy quantities): Array of the total mass enclosed for each corresponding distance in 
            Radius array. Astropy units of Msun. 
        """
        
        HaloMasses = self.MassEnclosed(Radius, 1) # compute halo mass enclosed by each distance in Radius
        DiskMasses = self.MassEnclosed(Radius, 2) # compute disk mass enclosed by each distance in Radius
        
        if self.gname != 'M33': # note: M33 has no bulge 
            BulgeMasses = self.MassEnclosed(Radius, 3) # compute bulge mass enclosed by each distance in Radius
        else:
            BulgeMasses = np.zeros(len(Radius)) * u.Msun # construct empty array for M33 
            
        TotalMasses = HaloMasses + DiskMasses + BulgeMasses # add up the mass arrays to get total mass at each distance
        
        # return array of total masses in units of Msun
        return TotalMasses
    
    def HernquistMass(self, Radius, a, m_halo): 
        """ 
        Method that defines the Hernquist 1990 mass profile 
        
        INPUTS:
            Radius: numpy array of astropy quantities
                Array of distances from the COM position in astropy units of kpc.
            a: astropy quantity
                scale radius of the Hernquist profile in kpc
            m_halo: astropy quantity 
                total halo mass in units of Msun 

        OUTPUTS:
            mass:  numpy array of astropy quantities 
                total mass within each distance in Radius in Msun
        """
       
        # Compute mass according to Hernquist profile equation 
        # NOTE: mass will be a numpy array because Radius is a numpy array
        mass = (m_halo * Radius**2) / (a + Radius)**2 

        # return Hernquist mass profile as array in units Msun
        return mass
    
    def CircularVelocity(self, Radius, ptype):
        """
        Method that computes the circular velocity of an object due to the mass enclosed in a galaxy within a certain radius.
        This velocity is ONLY due to the enclosed mass of a specific particle type. 
        
        INPUTS:
            Radius (numpy array of astropy quantities): Array of distances from the COM position in astropy units of kpc.
            ptype (int): Integer corresponding to a particle type (1 = Halo, 2 = Disk, 3 = Bulge)
            
        OUTPUTS:
            VCircs (numpy array of astropy quantities): Array of circular velocities at each distance in Radius.
                Units of km/s. 
        """
        
        # Get array of the masses enclosed within each distance in Radius for the given ptype 
        masses = self.MassEnclosed(Radius, ptype)
        
        # Construct empty numpy array to store circular velocities 
        VCircs = np.zeros(len(Radius))
        
        # Loop over each distance in Radius 
        for i in range(len(Radius)): 
            m = masses[i] # mass enclosed in Msun
            r = Radius[i] # distance in kpc 
            
            # Calculate circular velocity in km/s
            v_circ = np.sqrt(G*m/r)
            
            # Assign to corresponding index in VCircs 
            VCircs[i] = v_circ.value # note: v_circ is an astropy quantity and must be stored as a float first
             
        # Return array of the circular velocities rounded to 2 decimal places 
        return np.round(VCircs, 2) * u.km/u.s # assign units of km/s back 
     
    def TotalCircularVelocity(self, Radius):
        """
        Method that computes the circular velocity of an object due to the total mass enclosed in a galaxy 
        within a certain radius.
        
        INPUTS:
            Radius (numpy array of astropy quantities): Array of distances from the COM position in astropy units of kpc.
            
        OUTPUTS:
            TotalVCircs (numpy array of astropy quantities): Array of circular velocities at each distance in Radius.
                Units of km/s. 
        """
        
        # Get array of the total mass enclosed within each distance in Radius
        TotalMasses = self.TotalMassEnclosed(Radius)
        
        # Construct empty numpy array to store circular velocities 
        TotalVCircs = np.zeros(len(Radius))
        
        # Loop over each distance in Radius 
        for i in range(len(Radius)): 
            m = TotalMasses[i] # total mass enclosed in Msun 
            r = Radius[i] # radial distance in kpc 
            
            # Calculate circular velocity in km/s
            total_v_circ = np.sqrt(G*m/r)
            
            # Assign to corresponding index in TotalVCircs
            TotalVCircs[i] = total_v_circ.value # note: total_v_circ is an astropy quantity and must be stored as a float first
        
        # Return array of circular velocities in km/s
        return TotalVCircs * u.km / u.s # reassign units
    
    
    
    def HernquistVCirc(self, Radius, a, m_halo):
        """
        Method that computes a circular velocity due to the Hernquist mass enclosed in a galaxy within a given radius. 
        
        INPUTS:
            Radius: numpy array of astropy quantities
                Array of distances from the COM position in astropy units of kpc.
            a: astropy quantity
                scale radius of the Hernquist profile in kpc
            m_halo: astropy quantity 
                total halo mass in units of Msun 

        OUTPUTS:
            VCircs:  numpy array of astropy quantities 
                Circular velocities at each distance in Radius according to Hernquist profile. Units of km/s.
            
        """
        
        # Get array of the Hernquist masses at each distance in Radius. Units of Msun
        HernquistMasses = self.HernquistMass(Radius, a, m_halo)
        
        # Calculate circular velocity at each distance due to Hernquist mass enclosed. Units of km/s 
        # NOTE: Both HernquistMasses and Radius are arrays so VCircs will be an array 
        VCircs = np.sqrt(G * HernquistMasses / Radius)
        
        # Return array of circular velocities (km/s) rounded to 2 decimal places 
        return np.round(VCircs, 2)