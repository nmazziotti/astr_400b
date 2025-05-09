# Nicolas Mazziotti

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import os

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules from previous work 
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from CosmologicalTools import CosmologicalTools

# for contours
import scipy.optimize as so 

# for ellipse fitting 
from photutils.isophote import Ellipse, EllipseGeometry

# Define the benchmark cosmology at z =0
# Planck 2016 results. XIII. Cosmological parameters   
# Planck Collaboration+2016 A&A 594 13  Table 4, column 2 

OmegaM0_planck = 0.308   # Matter Density Parameter
OmegaR0_planck = 8.24e-5  # Radiation Density Parameter
OmegaL0_planck = 0.692  # Dark Energy Density Parameter
h_planck = 0.6781   # Hubble Constant  100 h km/s/Mpc

# Initialize Instances of the Cosmological Tools Class. 
BenchMark = CosmologicalTools(OmegaM0_planck,OmegaR0_planck,OmegaL0_planck,h_planck)

# Set Hubble constant to be value at present day (z=0)
H0 = BenchMark.HubbleParameter(0)

# Convert G to Mpc^3 /( Msun s^2 )
G = G.to(u.Mpc**3/u.Msun/u.s**2)

# Define critical density of the univerese in Msun/kpc^3
# Will use later on to calculate R200
rho_crit = (3*H0**2/(8*np.pi*G)).to(u.Msun/u.kpc**3)


# Taken from Lab 7
def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) 
            such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) 
            such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    
    # normalize the angular momentum vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to
    # z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel

# My class that I created
class HaloShape:
    '''
    Class that creates projections of M33's dark matter halo particle density and fits ellipses to the measure the halo shape
    '''

    def __init__(self, galaxy, snap):
        """
        PARAMETERS:
            - galaxy (string): Name of galaxy (e.g. 'M33')
            - snap (int): Desired simulation snapshot from 0 to 799
        """

        # Directory on nimoy where the high resolution data files are stored
        HighResDir = '../../../../astr400b/HighRes/M33/'
        
        # Determine Filename
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # create filenames
        self.filename = '%s_'%(galaxy) + ilbl + '.txt' # filename without directories 
        self.filesource = HighResDir + self.filename # filename with directories

        self.galaxy = galaxy # Store galaxy name 
        self.snap = snap  # Store snapshot 

        # Create a COM of object for M33 Halo (particle type=1) Using Code from Homework 4

        COMD = CenterOfMass(self.filesource,1)

        self.time = COMD.time
        # Compute COM of M33 using halo particles
        COMP = COMD.COM_P(0.1)
        COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])

        # Determine positions of halo particles relative to COM 
        xD = COMD.x - COMP[0].value 
        yD = COMD.y - COMP[1].value 
        zD = COMD.z - COMP[2].value 
        
        # total magnitude
        rtot = np.sqrt(xD**2 + yD**2 + zD**2)
        
        # Determine velocities of halo particles relative to COM motion
        vxD = COMD.vx - COMV[0].value 
        vyD = COMD.vy - COMV[1].value 
        vzD = COMD.vz - COMV[2].value 
        
        # total velocity 
        vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)
        
        # Arrays for r and v 
        r = np.array([xD,yD,zD]).T # transposed 
        v = np.array([vxD,vyD,vzD]).T

        # compute the rotated position and velocity vectors and store as global variables 
        self.rn, self.vn = RotateFrame(r, v)

    def ComputeR200(self, max_r=200):
        """
        Method to compute the R200 distance in kpc for the HaloShape object

        INPUTS:
            - max_r (int/float): Max radial separation in kpc from galaxy center of mass. Default is 200 kpc.
        """

        # Radial distance range in kpc 
        r = np.arange(0.25, max_r, 0.25) * u.kpc

        # Create MassProfile object for this galaxy
        Gal = MassProfile(self.galaxy, self.snap)

        # Compute DM mass enclosed at each radial distance 
        HaloMass = Gal.massEnclosed(1, r.value)

        # Convert masses into average densities by dividing each mass by a sphere 
        densities = 3*HaloMass/ (4*np.pi * r**3)

        # Find absolute difference between each density and 200 times the critical density 
        diff = abs(densities - 200 * rho_crit)

        # Find index where the minimum difference is 
        index = np.where(diff == np.min(diff))[0][0]

        # Set R200 to be the radial distance at this index in kpc
        r200 = r[index]

        return r200

    def GenerateProjection(self, axis1, axis2):
        """
        Method that generates a projection of the DM halo particle density in the rotated coordinate frame. Axes are 
        labeled as 0 for x, 1 for y, and 2 for z.

        INPUTS:
            - axis1 (int): Axis to display on the vertical axis 
            - axis2 (int): Axis to display on the horizontal axis 

        OUTPUTS:
            - h (2d numpy array): Projected histogram of particle density 
            - r200 (astropy quantity): R200 for this galaxy in kpc
        """

        fig, ax= plt.subplots(figsize=(12, 10))

        r200 = self.ComputeR200().value # Compute r200 at this snap and use to define box to make projection 

        # Create 2d histogram and store image data to variable called h 
        h = plt.hist2d(self.rn[:, axis1], self.rn[:, axis2], range = ((-r200, r200), (-r200, r200)), bins=75, norm=LogNorm(), cmap='cividis')
        
        # Hide figure
        plt.close(fig)

        # Return histogram data, r200 to use later in ellipse fitting 
        return h, r200

    def PlotEllipses(self, axis1, axis2, sma=10):
        """
        Method to plot fitted isodensity ellipses over 2D histogram projection of halo particle density 
        
        Inputs:
            - img_data (2D numpy array): Output 2D array from plt.hist2d 
            - sma (int, float): Initial guess of the semimajor axis length (in pixels) of the core region.
                Default value of 10 pixels.
            
        Outputs:
            Plots ellipses on top of projected image (axis labels not included)
        """

        # First generate projection image 
        hist, r200 = self.GenerateProjection(axis1, axis2)
        img_data = hist[0] # image is hist[0], horizontal coordinate information is hist[1], vertical coordinate information is hist[2] 

        # pixel bin size in kpc 
        x_binsize = (hist[2][1] - hist[2][0]) # same pixel scale for all 3 projections
        y_binsize = (hist[1][1] - hist[1][0])

        # Need to log10 the density data so that photutils can create ellipses based on a log scale 
        # Will not work without taking the log of the density data 
        fits_data = np.log10(img_data, where=(img_data > 0) ) # ignore where bin density is 0
        # fits_data will be used to perform the ellipse fitting 
        
        # Find where the highest density bin is and treat that as the object center (in pixel coordinates)
        x0 = np.where(img_data == np.nanmax(img_data))[1][0]
        y0 = np.where(img_data == np.nanmax(img_data))[0][0]
    
        
        ny, nx = fits_data.shape # Image size 
        
        # Initial guess of the ellipse geometry of the core region (sma will vary based on binning)
        geometry = EllipseGeometry(x0=x0, y0=y0, sma=sma, eps=0.3, pa=0.0, )
        
        # Create photutils Ellipse object using the fits_data and ellipse geometry 
        ellipse = Ellipse(fits_data, geometry)
        
        # Use photutils method fit_image to fit ellipses to the data out to a maximum sma length of 50 pixels 
        isolist = ellipse.fit_image(maxsma=40)
    
        #Useful plot that displays sma versus ellipticity to quantify the shape 
        fig, ax = plt.subplots(ncols=1, figsize=(10,5))
        
        # Multiply by pixel binszie to get semimajor axis in kpc 
        ax.errorbar(isolist.sma*y_binsize, isolist.eps, yerr=isolist.ellip_err, fmt='o', markersize=4)
        ax.set_xlabel('Semimajor axis length [kpc]')
        ax.set_ylabel('Ellipticity')
        ax.set_ylim(0,1) # 0 <= ellipticity <= 1
        ax.axvline(r200, c='magenta', linestyle=':', label=r"M33 $R_{200}$") # Plot dotted line at R200
        ax.axvline(r200/2, c='magenta', linestyle='--', label=r"M33 $R_{200}$/2") # Plot dashed line at R200/2
        ax.axvline(r200/4, c='magenta', linestyle='-.', label=r"M33 $R_{200}$/4") # Plot dotted-dashed line at R200/4
        ax.legend()

        # Set title to be snapshot time for analyzing time evolution 
        plt.title(np.round(self.time, 3))

        # Create new directory (if does not exist) to save ellipse plot to 
        save_ellipses_dir = f'./plots/{axis1}{axis2}/ellipses/' 
        if not os.path.isdir(save_ellipses_dir):
            os.makedirs(save_ellipses_dir)

        # filename will be M33_$$$.png, where $$$ is the snapshot number 
        fig.savefig(save_ellipses_dir + self.filename[:-4] + '.png', bbox_inches='tight')

        # Prevent figure from being plotted in jupyter notebook 
        plt.close(fig)

        # Create new figure to plot density projection with ellipse contours 
        fig2, ax2 = plt.subplots(figsize=(12, 10))

        # if the projection data doesn't work, the maximum density value will be 0. this usually occurs for the later snapshots    
        if np.max(img_data) == 0:
            plt.close(fig2) # close figure 
            raise ValueError # manually raise ValueError so that this will be considered a failed snapshot 
        else:
            # if projection data is successful, proceed 

            # show img_data in log scale
            sm = ax2.imshow(img_data, norm=LogNorm(),  cmap='cividis', interpolation='none')
            
            # Plot all of the fitted ellipses 
            smas = np.arange(0, max(isolist.sma), 1) # define array of semimajor axis lengths (in pixels) produced by the fitting process
            isos = np.empty(len(smas), dtype=object) # store isodensity contours to this array 
            for i in range(len(smas)):
                sma = smas[i] # current semimajor axis 
                iso = isolist.get_closest(sma) # get isodensity contour closest to this semimajor axis value 
                isos[i] = iso
                x, y, = iso.sampled_coordinates() # get x and y coordinates for this ellipse 
                plt.plot(x, y, color='lime', linewidth=0.5, linestyle='-') # plot on projection 

            # Initialize empty arrays to store ellipticity, position angle, and semimajor axis values at R200, R200/2, and R200/4
            r200_eps = np.zeros(3)
            r200_pa = np.zeros(3)
            r200_sma = np.zeros(3)

            # Get isodensity contours at the semimajor axis distances closest to R200, R200/2, and R200/4 in pixels 
            r200_smas = np.array([r200/y_binsize, r200/2/y_binsize, r200/4/y_binsize])
            for i in range(len(r200_smas)):
                sma = r200_smas[i]
                iso = isolist.get_closest(sma)
                r200_eps[i] = iso.eps # store ellipticity 
                r200_pa[i] = iso.pa # store position angle 
                r200_sma[i] = iso.sma * y_binsize # store semimajor axis in kpc 

                # plot these contours again on the projection but make them stand out with a magenta color 
                x, y, = iso.sampled_coordinates() 
                plt.plot(x, y, color='magenta', linewidth=1, linestyle='-')

            # conventions for axis labels according to axis1 and axis2 pairs 
            if axis1 == 0 and axis2 == 1:
                plt.ylabel('x [pixels]')
                plt.xlabel('y [pixels]')
            elif axis1 == 1 and axis2 == 2:
                plt.ylabel('y [pixels]')
                plt.xlabel('z [pixels]')
            elif axis1 == 0 and axis2 == 2:
                plt.ylabel('x [pixels]')
                plt.xlabel('z [pixels]')

            # Set title to be snapshot time for analyzing time evolution 
            plt.title(np.round(self.time, 3))

            # Add colorbar to plot (will be in log scale)
            cbar = plt.colorbar(sm)
            cbar.set_label("Number of DM particle per bin", fontsize=15)

            # Save to directory specifically for this projection in folder called projections 
            # if doesn't exist, create folders 
            save_projections_dir = f'./plots/{axis1}{axis2}/projections/' 
            if not os.path.isdir(save_projections_dir):
                os.makedirs(save_projections_dir)

            # filename will be M33_$$$.png, where $$$ is the snapshot number 
            fig2.savefig(save_projections_dir + self.filename[:-4] + '.png', bbox_inches='tight')

            # Prevent figure from being plotted in jupyter notebook 
            plt.close(fig2)

            # Also want to know semiminor axis at R200, R200/2, and R200/4, so calculate using ellipticity and semimajor axis in kpc 
            r200_smb = r200_sma * np.sqrt(1 - r200_eps)

            # Return isodensity contour objects and arrays of ellipse information
            return isolist, r200_eps, r200_pa, r200_sma, r200_smb