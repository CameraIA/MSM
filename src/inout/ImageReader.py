"""
.. module:: ImageReader
   :synopsis: This module is responsible for reading images as sequences of slices or unique stacks

.. moduleauthor:: T. Perciano

"""

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

from scipy import ndimage
from scipy import misc
from skimage import io
from colorama import init
from termcolor import colored
from ..util.util import create_dir

import os, re
import numpy as np
import shutil


from ..util.util import listdir_fullpath, create_dir

class ImageReader:
    """ An ImageReader

    :param input_dir:   Name of input directory
    :type input_dir:    str
    :param input_type:  Type of input
    :type input_type:   int
    :param fisrt_slice: First slice to be processed
    :type first_slice:  int
    :param last_slice: Last slice to be processed
    :type last_slice:  int
    :param create_dat:  If a .dat version of the image files should be created
    :type create_dat:   bool

    """
    def __init__(self, input_dir,input_type,first_slice,last_slice,create_dat,in_memory,run_number):
        """ ImageReader constructor

        """
        if input_type==0:
            self.input_dir = input_dir
        else:
            input_dir_split = input_dir.split("/")
            self.input_dir = "/".join(input_dir_split[:-1])+"/"
            self.input_file = input_dir
        self.create_dat = create_dat
        self.in_memory = in_memory
        self.run_number = run_number
        if self.create_dat:
            if input_type==0:
                self.input_dat_dir = os.path.join(self.input_dir,"msmcam_run_"+str(self.run_number), 'dat')
                create_dir(self.input_dat_dir)
            else:
                self.input_dat_dir = os.path.join(self.input_dir,"msmcam_run_"+str(self.run_number), 'dat')
                create_dir(self.input_dat_dir)

        if input_type==0:
            self.filenames = sorted(listdir_fullpath(str(self.input_dir)))
            """ List of file names for the input image """
            if '.DS_Store' in self.filenames[0]:
                self.filenames = self.filenames[1:]
            self.filenames = self.filenames[first_slice:last_slice]

        self.first_slice = first_slice
        self.last_slice = last_slice
        self.num_slices = self.last_slice-self.first_slice
        self.image_type = input_type
        self.image = None
        """ Input image """

    def getImage(self):
        """
        Method to return the input image

        :returns    Input image
        :rtype:     ndimage
        """
        return self.image

    def getImageFilenames(self):
        """
        Method to return the input image filenames

        :returns    Input image filenames
        :rtype:     list
        """
        return self.filenames

    def getNumSlices(self):
        """
        Method to return the number of slices

        :returns:   Number of slices of input image
        :rtype:     int
        """
        return self.num_slices

    def read(self):
        """
        Main procedure to read images
        """
        first_image = None
        if self.image_type==0:
            fname = self.filenames[0]
            try:
                first_image = io.imread(str(fname))
            except:
                print("Error reading image")
        if self.image_type==1:
            fname = self.input_file
            first_image = io.imread(str(fname))

        #: If image type is a sequence
        if self.image_type==0:
            y, x = first_image.shape
            z = self.num_slices
            if self.in_memory:
                temp_image = np.zeros((z,y,x), dtype=np.uint8)
            j = self.first_slice
            for i in range(z):
                img_unknown_type = io.imread(str(self.filenames[i]))
                if self.in_memory:
                    temp_image[i,:,:] = (255.*img_unknown_type/np.max(img_unknown_type)).astype(np.uint8)
                    misc.imsave(str(self.filenames[i]),temp_image[i,:,:])
                else:
                    temp_image = (255.*img_unknown_type/np.max(img_unknown_type)).astype(np.uint8)
                    misc.imsave(str(self.filenames[i]),temp_image)
                if self.create_dat:
                    name = str(self.input_dat_dir)+"/Slice"+"{0:0>4}".format(j)+".dat"
                    with open(name, 'wb') as f:
                        if self.in_memory:
                            temp_image[i,:,:].astype(np.uint8).tofile(f)
                        else:
                            temp_image.astype(np.uint8).tofile(f)
                print(colored('Read image '+str(j),'yellow','on_grey', attrs=['bold']))
                j = j+1

            if self.in_memory:
                self.image = temp_image

        #: If image type is a stack
        if self.image_type==1:
            new_list = list()
            first_image = first_image[self.first_slice:self.last_slice,:,:]
            z, y, x = first_image.shape
            self.image = (255.*first_image/np.max(first_image)).astype(np.uint8)
            j = self.first_slice
            for i in range(z):
                name = str(self.input_dir)+"/Slice"+"{0:0>4}".format(j)+".tif"
                misc.imsave(name,self.image[i,:,:])
                new_list.append(name)
                if self.create_dat:
                    name = str(self.input_dat_dir)+"/Slice"+"{0:0>4}".format(j)+".dat"
                    with open(name, 'wb') as f:
                        self.image[i,:,:].astype(np.uint8).tofile(f)
                j = j+1

            if self.in_memory is False:
                self.image = None
                
            self.filenames = new_list
