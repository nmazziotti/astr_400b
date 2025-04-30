# Nicolas Mazziotti

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile

# for contours
import scipy.optimize as so 

# for ellipse fitting 
from photutils.isophote import Ellipse, EllipseGeometry


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

    def ComputeR500(self, max_r=100.5):
        """
        Method to compute the R500 distance in kpc for the HaloShape object

        INPUTS:
            - max_r (int, float): Max radial separation in kpc from galaxy center of mass. Default is 100.5 kpc.
        """

        # Radial distance range in kpc 
        r = np.arange(0.25, max_r, 0.25) * u.kpc

        # Create MassProfile object for this galaxy
        Gal = MassProfile(self.galaxy, self.snap)

        # Compute DM mass enclosed at each radial distance 
        HaloMass = Gal.massEnclosed(1, r.value)

        # Convert masses into average densities by dividing each mass by a sphere 
        densities = 3*HaloMass/ (4*np.pi * r**3)

        # Define critical density of the univerese in Msun/kpc^3
        rho_crit = (9.47e-27 * u.kg/u.m**3).to(u.Msun/u.kpc**3)

        # Find absolute difference between each density and 500 times the critical density 
        diff = abs(densities - 500 * rho_crit)

        # Find index where the minimum difference is 
        index = np.where(diff == np.min(diff))[0][0]

        # Set R500 to be the radial distance at this index in kpc
        r500 = r[index]

        return r500

    def GenerateProjection(self, axis1, axis2):
        """
        Method that generates a projection of the DM halo particle density in the rotated coordinate frame. Axes are 
        labeled as 0 for x, 1 for y, and 2 for z.

        INPUTS:
            - axis1 (int): Axis to display on the vertical axis 
            - axis2 (int): Axis to display on the horizontal axis 

        OUTPUTS:
            - h (2d numpy array): Projected histogram of particle density 
            - r500 (astropy quantity): R500 for this galaxy in kpc
        """

        fig, ax= plt.subplots(figsize=(12, 10))

        r500 = self.ComputeR500().value # Compute r500 at this snap and use to define box to make projection 
        h = plt.hist2d(self.rn[:, axis1], self.rn[:, axis2], range = ((-r500, r500), (-r500, r500)), bins=75, norm=LogNorm(), cmap='cividis')
        
        cbar = plt.colorbar()
        cbar.set_label("Number of DM particle per bin", fontsize=15)
        
        
        #adjust tick label font size
        label_size = 10
        matplotlib.rcParams['xtick.labelsize'] = label_size 
        matplotlib.rcParams['ytick.labelsize'] = label_size
        
        # Hide figure 
        plt.close(fig)

        return h, r500

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
        
        # Need to log10 the density data so that photutils can create ellipses based on a log scale 
        # Will not work without taking the log of the density data 

        hist, r500 = self.GenerateProjection(axis1, axis2)
        img_data = hist[0]

        x_binsize = (hist[2][1] - hist[2][0]) # same for all 3 projections
        y_binsize = (hist[1][1] - hist[1][0])

        fits_data = np.log10(img_data, where=(img_data != 0) ) # ignore where bin density is 0
        
        # Find where the highest density bin is and treat that as the object center (in pixel coordinates)
        x0 = np.where(img_data == np.max(img_data))[1][0]
        y0 = np.where(img_data == np.max(img_data))[0][0]
    
        
        ny, nx = fits_data.shape # Image size 
        
        # Initial guess of the ellipse geometry of the core region (sma will vary based on binning)
        geometry = EllipseGeometry(x0=x0, y0=y0, sma=sma, eps=0.3, pa=0.0, )
        
        # Create photutils Ellipse object using the fits_data and ellipse geometry 
        ellipse = Ellipse(fits_data, geometry)
        
        # Use photutils method fit_image to fit ellipses to the data out to a maximum sma length of 50 pixels 
        isolist = ellipse.fit_image(maxsma=40)
    
        #Useful plot that displays sma versus ellipticity and position angle to quantify the shape 
        fig, ax = plt.subplots(ncols=2, figsize=(10,5))
        ax[0].errorbar(isolist.sma*y_binsize, isolist.eps, yerr=isolist.ellip_err, fmt='o', markersize=4)
        ax[0].set_xlabel('Semimajor axis length [kpc]')
        ax[0].set_ylabel('Ellipticity')
        ax[0].set_ylim(0,1)
        ax[0].axvline(r500, c='magenta', linestyle=':', label=r"M33 $R_{500}$")
        ax[0].axvline(r500/2, c='magenta', linestyle='--', label=r"M33 $R_{500}$/2")
        ax[0].axvline(r500/4, c='magenta', linestyle='-.', label=r"M33 $R_{500}$/4")
        ax[0].legend()
        
        ax[1].errorbar(isolist.sma*y_binsize, isolist.pa, yerr=isolist.pa_err, fmt='o', markersize=4)
        ax[1].set_xlabel('Semimajor axis length [kpc]')
        ax[1].set_ylabel('Position Angle')
        ax[1].axvline(r500, c='magenta', linestyle=':', label=r"M33 $R_{500}$")
        ax[1].axvline(r500/2, c='magenta', linestyle='--', label=r"M33 $R_{500}$/2")
        ax[1].axvline(r500/4, c='magenta', linestyle='-.', label=r"M33 $R_{500}$/4")
        ax[1].legend()
        
        #plt.errorbar(isolist.sma, isolist.pa, yerr=isolist.pa_err, fmt='o', markersize=4)

        plt.title(np.round(self.time, 3))
        plt.savefig(f'./plots/{axis1}{axis2}/ellipses/' + self.filename[:-4] + '.png', bbox_inches='tight')
        plt.close(fig)
        
    
        fig, ax = plt.subplots(figsize=(12, 10))
        sm = ax.imshow(img_data, norm=LogNorm(),  cmap='cividis', interpolation='none')
        w = nx * 0.05 # Show image box centered on the halo and is 10% the total image size in pixel length 
    #     ax.set_xlim(x0 - w, x0 + w)
    #     ax.set_ylim(y0 - w, y0 + w)
        
        # Plot the fitted ellipses 
        isos = []
        smas = np.arange(0, max(isolist.sma), 1)
        for sma in smas:
            iso = isolist.get_closest(sma)
            isos.append(iso)
            x, y, = iso.sampled_coordinates()
            plt.plot(x, y, color='lime', linewidth=0.5, linestyle='-')

        r500_eps = []
        r500_pa = []
        r500_sma = []
        
        for sma in [r500/y_binsize, r500/2/y_binsize, r500/4/y_binsize]:
            iso = isolist.get_closest(sma)
            r500_eps.append(iso.eps)
            r500_pa.append(iso.pa)
            r500_sma.append(iso.sma * y_binsize)
            
            x, y, = iso.sampled_coordinates()
            plt.plot(x, y, color='magenta', linewidth=1, linestyle='-')

        plt.title(np.round(self.time, 3))
        cbar = plt.colorbar(sm)
        cbar.set_label("Number of DM particle per bin", fontsize=15)
        
        plt.savefig(f'./plots/{axis1}{axis2}/projections/' + self.filename[:-4] + '.png', bbox_inches='tight')
        plt.close(fig)

        r500_eps = np.array(r500_eps)
        r500_pa = np.array(r500_pa)
        r500_sma = np.array(r500_sma)
        r500_smb = r500_sma * np.sqrt(1 - r500_eps)
        
        return isolist, r500_eps, r500_pa, r500_sma, r500_smb
            