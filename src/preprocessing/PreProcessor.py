"""
.. module:: Preprocessor
   :synopsis: This module is responsible for filtering input data

.. moduleauthor:: T. Perciano

"""

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage.transform import downscale_local_mean
from skimage import restoration as rt
from skimage import io, img_as_float, img_as_ubyte, segmentation, filters
from scipy import misc
import os
import shutil
import numpy as np
from colorama import init
from termcolor import colored

from ..util.util import create_dir, listdir_fullpath

class PreProcessor:
    """ A PreProcessor

    :param input_image:         Image to be filtered
    :type input_image:          ndimage or list (names of files)
    :param input_settings:		Dictionary with input settings
    :type input_settings:		dic
    :param preproc_settings:	Dictionary with preprocessing settings
    :type preproc_settings:		dic

    """

    def __init__(self, input_image, input_settings, preproc_settings, run_number):
        """ Initialize Preprocessor thread

        """
        self.input_image = input_image
        """ Image to be filtered """
        self.input_type = input_settings["InputType"]
        if self.input_type==0:
            self.input_dir = input_settings["InputDir"]
        else:
            input_dir = input_settings["InputDir"]
            input_dir_split = input_dir.split("/")
            self.input_dir = "/".join(input_dir_split[:-1])+"/"
            self.input_file = input_dir
        self.down_scale = preproc_settings["DownsampleScale"]
        """ Flag indicating if image should be downsampled """
        self.filtered = self.input_image
        self.filtered_filenames = list()
        """ Filtered data """
        self.usemedian = preproc_settings["UseMedian"]
        """ Flag indicating if median filter should be used """
        self.usebilateral = preproc_settings["UseBilateral"]
        """ Flag indicating if bilateral filter should be used """
        self.medianmasksize = preproc_settings["MedianMaskSize"]
        """ Size of median filter mask """
        self.sigmaspatial = preproc_settings["BilateralSigmaSpatial"]
        """ Sigma spatial parameter for bilateral filter """
        self.sigmacolor = preproc_settings["BilateralSigmaColor"]
        """ Sigma color parameter for bilateral filter """
        self.in_memory = input_settings["InMemory"]
        self.run_number = run_number
        if self.in_memory is False:
            self.dat_dir = self.input_dir+ "msmcam_run_"+str(self.run_number)+'/preproc/dat'
            """ Output directory for preprocessed data in .dat format """
            self.tif_dir = self.input_dir+ "msmcam_run_"+str(self.run_number)+'/preproc/tif'
            """ Output directory for preprocessed data in .tif format """
            create_dir(self.dat_dir)
            create_dir(self.tif_dir)

        self.first_slice = int(input_settings["FirstSlice"])
        self.last_slice = int(input_settings["LastSlice"])
        self.num_slices = self.last_slice-self.first_slice
         
    def getFiltered(self):
        """ Returns filtered image

        :returns:   Filtered image
        :rtype:     ndimage

        """
        return self.filtered, self.filtered_filenames

    def downsample(self):
        """ Downsamples image

        :returns:   Downsampled image
        :rtype:     ndimage

        """

        down = downscale_local_mean(self.filtered, (self.down_scale,self.down_scale,self.down_scale))
        down = down.astype(np.uint8)

        return(down)

    def median(self,image):
        """ Applied median filter to image
        
        """
        return(ndi.filters.median_filter(image,size=self.medianmasksize))

    def bilateral(self,image):
        """ Applies bilateral filter to image

        """
        return(img_as_ubyte(rt.denoise_bilateral(image,sigma_color=self.sigmacolor,sigma_spatial=self.sigmaspatial,multichannel=False)))

    def nlfilter(self):
        """ Non-linear filtering """
        #TODO

    def save(self,image,nslice):
        """ Saves filtered image to disk """
        name = str(self.tif_dir)+"/Slice"+"{0:0>4}".format(nslice)+".tif"
        misc.imsave(name,image)
        self.filtered_filenames.append(name)

        name = str(self.dat_dir)+"/Slice"+"{0:0>4}".format(nslice)+".dat"
        with open(name, 'wb') as f:
            temp = image.astype(np.uint8)
            temp.tofile(f)

    def process(self):
        """ Complete filtering pipeline """
        print(colored('Denoising images...','green','on_grey', attrs=['bold']))
        if self.in_memory:
            self.filtered = self.downsample()
            
        j = self.first_slice
        for i in range(self.num_slices):
            filtered = None
            if self.usemedian:
                if self.in_memory:
                    self.filtered[i,:,:] = self.median(self.filtered[i,:,:])
                else:
                    image = io.imread(str(self.filtered[i]))
                    filtered = self.median(image)
            if self.usebilateral:
                if self.in_memory:
                    self.filtered[i,:,:] = self.bilateral(img_as_float(self.filtered[i,:,:]))
                else:
                    image = io.imread(str(self.filtered[i]))
                    if filtered is None:
                        filtered = self.bilateral(img_as_float(image))
                    else:
                        filtered = self.bilateral(img_as_float(filtered))
            if self.in_memory is False:
                if filtered is None:
                    filtered = io.imread(str(self.filtered[i]))
                self.save(filtered,j)
            print(colored('Finished image '+str(j),'yellow','on_grey', attrs=['bold']))
            j = j+1
            
