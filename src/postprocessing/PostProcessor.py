"""
.. module:: Postprocessor
   :synopsis: This module is responsible for applying postprocessing routines

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

from colorama import init
from termcolor import colored
from skimage import io, img_as_float, img_as_ubyte, img_as_bool, morphology
from scipy import misc
from ..util.util import create_dir, listdir_fullpath

class PostProcessor:
    """ A PostProcessor



    """

    def __init__(self, input_image, input_settings, postprocess_settings, run_number):
        """ Initialize Postprocessor thread

        """
        input_type = input_settings["InputType"]
        input_dir = input_settings["InputDir"]
        if input_type==1:
            input_dir_split = input_dir.split("/")
            input_dir = "/".join(input_dir_split[:-1])

        self.postprocessed = input_image
        self.postprocessed_filenames = list()
        self.in_memory = input_settings["InMemory"]
        self.run_number = run_number
        if self.in_memory is False:
            self.tif_dir = input_dir+ "/msmcam_run_"+str(self.run_number)+'/postproc/tif'
            """ Output directory for postprocessed data in .tif format """
            create_dir(self.tif_dir)

        self.first_slice = input_settings["FirstSlice"]
        self.last_slice = input_settings["LastSlice"]
        self.num_slices = self.last_slice-self.first_slice
        
        self.remove_small_objects = postprocess_settings["RemoveSmallObjects"]
        self.remove_small_holes = postprocess_settings["RemoveSmallHoles"]

    def save(self,image,nslice):
        """ Saves postprocessed image to disk """
        name = str(self.tif_dir)+"/Slice"+"{0:0>4}".format(nslice)+".tif"
        misc.imsave(name, image)
        self.postprocessed_filenames.append(name)
        
    def getResult(self):
        """ Returns postprocessed image

        :returns:   Postprocessed image
        :rtype:     ndimage

        """
        return(self.postprocessed)
        
    def getResultFilenames(self):
        """ Returns postprocessed image (filenames)

        :returns:   Postprocessed image (filenames)
        :rtype:     list

        """
        return(self.postprocessed_filenames)
       
    def process(self):
        """ Complete postprocessing pipeline """
        print(colored('Postprocessing images...','green','on_grey', attrs=['bold']))
          
        j = self.first_slice
        for i in range(self.num_slices):
            postprocessed = None
            if self.remove_small_objects:
                if self.in_memory:
                    self.postprocessed[i,:,:] = img_as_ubyte(morphology.remove_small_objects(img_as_bool(self.postprocessed[i,:,:]), minSize, connectivity=2))
                else:
                    image = io.imread(str(self.postprocessed[i]))
                    filtered = img_as_ubyte(morphology.remove_small_objects(img_as_bool(image), minSize, connectivity=2))
            if self.remove_small_holes:
                if self.in_memory:
                    self.postprocessed[i,:,:] = img_as_ubyte(morphology.remove_small_holes(img_as_bool(self.postprocessed[i,:,:]), min_size=100, connectivity=4, in_place=True))
                else:
                    image = io.imread(str(self.postprocessed[i]))
                    if postprocessed is None:
                        postprocessed = img_as_ubyte(morphology.remove_small_holes(img_as_bool(image), min_size=100, connectivity=4, in_place=True))
                    else:
                        postprocessed = img_as_ubyte(morphology.remove_small_holes(img_as_bool(postprocessed), min_size=100, connectivity=4, in_place=True))
            if self.in_memory is False:
                if postprocessed is None:
                    postprocessed = io.imread(str(self.postprocessed[i]))
                self.save(postprocessed,j)
            print(colored('Finished image '+str(j),'yellow','on_grey', attrs=['bold']))
            j = j+1
