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

    def __init__(self, galaxy, snap):
        HighResDir = '../../../../astr400b/HighRes/M33/'
        
        # Determine Filename
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # create filenames
        self.filename = '%s_'%(galaxy) + ilbl + '.txt'
        self.filesource = HighResDir + self.filename

        self.snap = snap 

        # Create a COM of object for M33 Halo (particle type=1) Using Code from Homework 4

        # Using HighRes version of M33_000.txt 
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

        # ADD HERE
        # compute the rotated position and velocity vectors
        self.rn, self.vn = RotateFrame(r, v)

    def GenerateProjection(self, axis1, axis2):
        # 1) Make plots 

        # M33 Disk Density 
        fig, ax= plt.subplots(figsize=(12, 10))
        
        # ADD HERE
        # plot the particle density for M33 using a 2D historgram
        # plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
        # cmap options: 
        # https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
        #   e.g. 'magma', 'viridis'
        # can modify bin number to make the plot smoother
        h = plt.hist2d(self.rn[:, axis1], self.rn[:, axis2], range = ((-50, 50), (-50, 50)), bins=75, norm=LogNorm(), cmap='cividis')
        
        cbar = plt.colorbar()
        cbar.set_label("Number of DM particle per bin", fontsize=15)
        
        
        #set axis limits
        # l = 0.5e3 # length of box to show in kpc
        # plt.ylim(-l,l)
        # plt.xlim(-l,l)
        
        #adjust tick label font size
        label_size = 10
        matplotlib.rcParams['xtick.labelsize'] = label_size 
        matplotlib.rcParams['ytick.labelsize'] = label_size
        
        
        # Save to a file
        #plt.savefig('Lab7_M31DMHalo.png', bbox_inches='tight')
        plt.close(fig)

        return h

    def PlotEllipses(self, axis1, axis2, sma=10):
        """
        Function to plot fitted isodensity ellipses over 2D histogram projection of halo particle density 
        
        Inputs:
            - img_data (2D numpy array): Output 2D array from plt.hist2d 
            - sma (int, float): Initial guess of the semimajor axis length (in pixels) of the core region.
                Default value of 10 pixels.
            
        Outputs:
            Plots ellipses on top of projected image (axis labels not included)
        """
        
        # Need to log10 the density data so that photutils can create ellipses based on a log scale 
        # Will not work without taking the log of the density data 

        hist = self.GenerateProjection(axis1, axis2)
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
        isolist = ellipse.fit_image(maxsma=37.5)
    
        #Useful plot that displays sma versus ellipticity and position angle to quantify the shape 
        fig, ax = plt.subplots(ncols=2, figsize=(10,5))
        ax[0].errorbar(isolist.sma*y_binsize, isolist.eps, yerr=isolist.ellip_err, fmt='o', markersize=4)
        ax[0].set_xlabel('Semimajor axis length [kpc]')
        ax[0].set_ylabel('Ellipticity')
        ax[0].set_ylim(0,1)
        ax[0].axvline(25, c='magenta', linestyle=':', label="M33 Scale Radius")
        ax[0].axvline(25/2, c='magenta', linestyle='--', label="Half of Scale Radius")
        ax[0].legend()
        
        ax[1].errorbar(isolist.sma*y_binsize, isolist.pa, yerr=isolist.pa_err, fmt='o', markersize=4)
        ax[1].set_xlabel('Semimajor axis length [kpc]')
        ax[1].set_ylabel('Position Angle')
        ax[1].axvline(25, c='magenta', linestyle=':', label="M33 Scale Radius")
        ax[1].axvline(25/2, c='magenta', linestyle='--', label="Half of Scale Radius")
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
            
        for sma in [25/y_binsize, 25/2/y_binsize]:
            iso = isolist.get_closest(sma)
            x, y, = iso.sampled_coordinates()
            plt.plot(x, y, color='magenta', linewidth=1, linestyle='-')

        plt.title(np.round(self.time, 3))
        cbar = plt.colorbar(sm)
        cbar.set_label("Number of DM particle per bin", fontsize=15)
        
        plt.savefig(f'./plots/{axis1}{axis2}/projections/' + self.filename[:-4] + '.png', bbox_inches='tight')
        plt.close(fig)
        return isolist
            