# This function has the sets the calibration and mask for the GIWAXS experiment of interest
# We find that 


import sys
import re # Regular expressions

import numpy as np
import matplotlib as mpl
from SciAnalysis.settings import * # we might not need this

mpl.rcParams['mathtext.fontset'] = 'cm'
import pylab as plt


# Calibration
################################################################################    
class Calibration(object):
    '''Stores aspects of the experimental setup; especially the calibration
    parameters for a particular detector. That is, the wavelength, detector
    distance, and pixel size that are needed to convert pixel (x,y) into
    reciprocal-space (q) value.
    
    This class may also store other information about the experimental setup
    (such as beam size and beam divergence).
    '''
    
    def __init__(self, wavelength_A=None, distance_m=None, pixel_size_um=None, incident_angle=0):
        
        self.wavelength_A = wavelength_A
        self.distance_m = distance_m
        self.pixel_size_um = pixel_size_um
        
        self.incident_angle = incident_angle
        self.sample_normal = None
        self._beam_positions = {}
        
        
        # Data structures will be generated as needed
        # (and preserved to speedup repeated calculations)
        self.clear_maps()
    
    
    # Experimental parameters
    ########################################
    
    def set_wavelength(self, wavelength_A):
        '''Set the experimental x-ray wavelength (in Angstroms).'''
        
        self.wavelength_A = wavelength_A
    
    
    def get_wavelength(self):
        '''Get the x-ray beam wavelength (in Angstroms) for this setup.'''
        
        return self.wavelength_A
    
        
    def set_energy(self, energy_keV):
        '''Set the experimental x-ray beam energy (in keV).'''
        
        energy_eV = energy_keV*1000.0
        energy_J = energy_eV/6.24150974e18
        
        h = 6.626068e-34 # m^2 kg / s
        c = 299792458 # m/s
        
        wavelength_m = (h*c)/energy_J
        self.wavelength_A = wavelength_m*1e+10
    
    
    def get_energy(self):
        '''Get the x-ray beam energy (in keV) for this setup.'''
        
        h = 6.626068e-34 # m^2 kg / s
        c = 299792458 # m/s
        
        wavelength_m = self.wavelength_A*1e-10 # m
        E = h*c/wavelength_m # Joules
        
        E *= 6.24150974e18 # electron volts
        
        E /= 1000.0 # keV
        
        return E
    
    
    def get_k(self):
        '''Get k = 2*pi/lambda for this setup, in units of inverse Angstroms.'''
        
        return 2.0*np.pi/self.wavelength_A
    
    
    def set_distance(self, distance_m):
        '''Sets the experimental detector distance (in meters).'''
        
        self.distance_m = distance_m
        
    
    def set_pixel_size(self, pixel_size_um=None, width_mm=None, num_pixels=None):
        '''Sets the pixel size (in microns) for the detector. Pixels are assumed
        to be square.'''
        
        if pixel_size_um is not None:
            self.pixel_size_um = pixel_size_um
            
        else:
            if num_pixels is None:
                num_pixels = self.width
            pixel_size_mm = width_mm*1./num_pixels
            self.pixel_size_um = pixel_size_mm*1000.0
        
        
    def set_beam_position(self, x0, y0, named=None):
        '''Sets the direct beam position in the detector images (in pixel 
        coordinates).'''
        
        if named is not None:
            self._beam_positions[named] = [x0, y0]
        else:
            self._beam_positions['default'] = [x0, y0]
            self.x0 = x0
            self.y0 = y0
        
        
    def use_beam_position(self, name):
        self.x0, self.y0 = self._beam_positions[name]
        
        
    def set_image_size(self, width, height=None):
        '''Sets the size of the detector image, in pixels.'''
        
        self.width = width
        if height is None:
            # Assume a square detector
            self.height = width
        else:
            self.height = height
    
    
    def get_q_per_pixel(self):
        '''Gets the delta-q associated with a single pixel. This is computed in
        the small-angle limit, so it should only be considered approximate.
        For instance, wide-angle detectors will have different delta-q across
        the detector face.'''
        
        if self.q_per_pixel is not None:
            return self.q_per_pixel
        
        c = (self.pixel_size_um/1e6)/self.distance_m
        twotheta = np.arctan(c) # radians
        
        self.q_per_pixel = 2.0*self.get_k()*np.sin(twotheta/2.0)
        
        return self.q_per_pixel
    
    
    def set_incident_angle(self, incident_angle=0, sample_normal=None):
        
        self.clear_maps() # Presumptively invalidate cached maps
        self.incident_angle = incident_angle
        if sample_normal is not None:
            self.sample_normal = sample_normal
    
    
    def set_angles(self, sample_normal=0, incident_angle=None):
        
        self.clear_maps() # Presumptively invalidate cached maps
        self.sample_normal = sample_normal
        if incident_angle is not None:
            self.incident_angle = incident_angle
    
    
    # Convenience methods
    ########################################
    def q_to_angle(self, q):
        '''Convert from q to angle (full scattering angle, 2theta, in degrees).'''
        kpre = 2.0*self.get_k()
        return np.degrees( 2.0*np.arcsin(q/kpre) )
    
    def angle_to_q(self, angle):
        '''Convert from scattering angle (full scattering angle, in degrees)
        to q-value (in inverse angstroms).'''
        kpre = 2.0*self.get_k()
        return kpre*np.sin(np.radians(angle/2))    
    
    def compute_qy(self, QX, QZ):
        '''Compute the (sometimes ignored!) qy component for the given (QX,QZ) matrices.'''
        
        k = self.get_k()
        
        
        alpha_f = np.arcsin(QZ/k)
        theta_f = np.arcsin( QX/(k*np.cos(alpha_f)) )
        
        QY = k*( np.cos(theta_f)*np.cos(alpha_f) - 1 )
        
        # Alternate computation:
        #QZk2 = (1 - np.square(QZ/k))
        #QY = k*( np.sqrt( 1 - np.square(QX/k)*(1/QZk2) )*np.sqrt(QZk2) - 1 )
        
        return QY
    
    # Maps
    ########################################
    
    def clear_maps(self):
        self.r_map_data = None
        self.q_per_pixel = None
        self.q_map_data = None
        self.angle_map_data = None
        
        self.qx_map_data = None
        self.qy_map_data = None
        self.qz_map_data = None
        self.qr_map_data = None

    
    def r_map(self):
        '''Returns a 2D map of the distance from the origin (in pixel units) for
        each pixel position in the detector image.'''
        
        if self.r_map_data is not None:
            return self.r_map_data

        x = np.arange(self.width) - self.x0
        y = np.arange(self.height) - self.y0
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        self.r_map_data = R
        
        return self.r_map_data
        
    
    def q_map(self):
        '''Returns a 2D map of the q-value associated with each pixel position
        in the detector image.'''

        if self.q_map_data is not None:
            return self.q_map_data
        
        c = (self.pixel_size_um/1e6)/self.distance_m
        twotheta = np.arctan(self.r_map()*c) # radians
        
        self.q_map_data = 2.0*self.get_k()*np.sin(twotheta/2.0)
        
        return self.q_map_data
        
    
    def angle_map(self):
        '''Returns a map of the angle for each pixel (w.r.t. origin).
        0 degrees is vertical, +90 degrees is right, -90 degrees is left.'''

        if self.angle_map_data is not None:
            return self.angle_map_data
        
        x = (np.arange(self.width) - self.x0)
        y = (np.arange(self.height) - self.y0)
        X,Y = np.meshgrid(x,y)
        #M = np.degrees(np.arctan2(Y, X))
        # Note intentional inversion of the usual (x,y) convention.
        # This is so that 0 degrees is vertical.
        #M = np.degrees(np.arctan2(X, Y))

        # TODO: Lookup some internal parameter to determine direction
        # of normal. (This is what should befine the angle convention.)
        M = np.degrees(np.arctan2(X, -Y))

        
        self.angle_map_data = M

        if self.sample_normal is not None:
            self.angle_map_data += self.sample_normal
        
        
        return self.angle_map_data
    
    
    def qx_map(self):
        if self.qx_map_data is not None:
            return self.qx_map_data
        
        self._generate_qxyz_maps()
        
        return self.qx_map_data    

    def qy_map(self):
        if self.qy_map_data is not None:
            return self.qy_map_data
        
        self._generate_qxyz_maps()
        
        return self.qy_map_data    

    def qz_map(self):
        if self.qz_map_data is not None:
            return self.qz_map_data
        
        self._generate_qxyz_maps()
        
        return self.qz_map_data    
    
    def qr_map(self):
        if self.qr_map_data is not None:
            return self.qr_map_data

        self._generate_qxyz_maps()

        return self.qr_map_data



    def _generate_qxyz_maps(self):
        
        # Conversion factor for pixel coordinates
        # (where sample-detector distance is set to d = 1)
        c = (self.pixel_size_um/1e6)/self.distance_m
        
        x = np.arange(self.width) - self.x0
        y = np.arange(self.height) - self.y0
        X, Y = np.meshgrid(x, y)
        R = np.sqrt(X**2 + Y**2)
        
        #twotheta = np.arctan(self.r_map()*c) # radians
        theta_f = np.arctan2( X*c, 1 ) # radians
        #alpha_f_prime = np.arctan2( Y*c, 1 ) # radians
        alpha_f = np.arctan2( Y*c*np.cos(theta_f), 1 ) # radians
        
        cos_inc = np.cos(np.radians(self.incident_angle))
        sin_inc = np.sin(np.radians(self.incident_angle))
        self.qx_map_data = self.get_k()*np.sin(theta_f)*np.cos(alpha_f)
        self.qy_map_data = self.get_k()*( np.cos(theta_f)*np.cos(alpha_f) - cos_inc ) # TODO: Check sign
        self.qz_map_data = -1.0*self.get_k()*( np.sin(alpha_f) + sin_inc ) 
        

        if self.sample_normal is not None:
            s = np.sin(np.radians(self.sample_normal))
            c = np.cos(np.radians(self.sample_normal))
            self.qx_map_data, self.qz_map_data = c*self.qx_map_data - s*self.qz_map_data, s*self.qx_map_data + c*self.qz_map_data
        
        self.qr_map_data = np.sign(self.qx_map_data)*np.sqrt(np.square(self.qx_map_data) + np.square(self.qy_map_data))


 # Mask
################################################################################    
class Mask(object):
    '''Stores the matrix of pixels to be excluded from further analysis.'''
    
    def __init__(self, infile=None, format='auto'):
        '''Creates a new mask object, storing a matrix of the pixels to be 
        excluded from further analysis.'''
        
        self.data = None
        
        if infile is not None:
            self.load(infile, format=format)
        
        
    def load(self, infile, format='auto', invert=False):
        '''Loads a mask from a a file. If this object already has some masking
        defined, then the new mask is 'added' to it. Thus, one can load multiple
        masks to exlude various pixels.'''
        
        if format=='png' or infile[-4:]=='.png':
            self.load_png(infile, invert=invert)
            
        elif format=='hdf5' or infile[-3:]=='.h5' or infile[-4:]=='.hd5':
            self.load_hdf5(infile, invert=invert)
            
        else:
            print("Couldn't identify mask format for %s."%(infile))
            
            
    def load_blank(self, width, height):
        '''Creates a null mask; i.e. one that doesn't exlude any pixels.'''
        
        # TODO: Confirm that this is the correct order for x and y.
        self.data = np.ones((height, width))
        
            
    def load_png(self, infile, threshold=127, invert=False):
        '''Load a mask from a PNG image file. High values (white) are included, 
        low values (black) are exluded.'''
        
        # Image should be black (0) for excluded pixels, white (255) for included pixels
        img = PIL.Image.open(infile).convert("L") # black-and-white
        img2 = img.point(lambda p: p > threshold and 255)
        data = np.asarray(img2)/255
        data = data.astype(int)
        
        if invert:
            data = -1*(data-1)
        
        if self.data is None:
            self.data = data
        else:
            self.data *= data
        
        
    def load_hdf5(self, infile, invert=False):
        
        with h5py.File(infile, 'r') as f:
            data = np.asarray( f['mask'] )

        if invert:
            data = -1*(data-1)
        
        if self.data is None:
            self.data = data
        else:
            self.data *= data

        
    def invert(self):
        '''Inverts the mask. Can be used if the mask file was written using the
        opposite convention.'''
        self.data = -1*(self.data-1)


    # End class Mask(object)
    ########################################