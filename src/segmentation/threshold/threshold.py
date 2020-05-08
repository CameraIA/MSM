"""
    .. module:: Threshold
    :synopsis: This module is responsible for image segmentation using any threshold algorithm
    
    .. moduleauthor:: T. Perciano
    
    """

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"


from ...util.util import create_dir, upscale, listdir_fullpath
from skimage import io, img_as_ubyte, morphology, measure
from scipy import misc
from colorama import init
from termcolor import colored
from ..SegmentationAlgorithm import SegmentationAlgorithm

import numpy as np
import skimage.filters as filt
import os

class Threshold(SegmentationAlgorithm):
    def __init__(self, image, input_settings, preproc_settings,seg_settings,run_number):
        super().__init__(image, input_settings, preproc_settings,seg_settings,run_number)

    def segment(self):
        """ Main segmentation process using a threshold algorithm 
            
        """
        
        input_type = self.input_settings["InputType"]
        input_dir = self.input_settings["InputDir"]
        if input_type==1:
            input_dir_split = input_dir.split("/")
            input_dir = "/".join(input_dir_split[:-1])
            input_file = input_dir
        
        masked = self.input_settings["Masked"]
        invert = self.seg_settings["Invert"]
        in_memory = self.input_settings["InMemory"]
        first_slice = self.input_settings["FirstSlice"]
        last_slice = self.input_settings["LastSlice"]
        alg = self.seg_settings["ThreshAlg"]
        num_slices = last_slice-first_slice
        if in_memory:
            self.segmented = np.zeros(self.image.shape,dtype=np.uint8)
        
        # create outout directories, if writing tifs
        if not in_memory:
            output_tif_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number),"res/thresh/tif")
            create_dir(output_tif_dir)

        #if in_memory:
        #    self.image = upscale(self.image, self.preproc_settings["DownsampleScale"])
        j = first_slice
        for i in range(num_slices):
            if in_memory:
                img = self.image[i,:,:]
            else:
                img = io.imread(str(self.image[i]))
            
            # Apply threshold
            if alg=='Otsu':
                img_res = img <= filt.threshold_otsu(img)
            elif alg=='Yen':
                img_res = img <= filt.threshold_yen(img)
            elif alg=='Isodata':
                img_res = img <= filt.threshold_isodata(img)
            elif alg=='Li':
                img_res = img <= filt.threshold_li(img)
            elif alg=='Minimum':
                img_res = img <= filt.threshold_minimum(img)
            elif alg=='Mean':
                img_res = img <= filt.threshold_mean(img)
            elif alg=='Niblack':
                img_res = img <= filt.threshold_niblack(img)
            elif alg=='Sauvola':
                img_res = img <= filt.threshold_sauvola(img)
            elif alg=='Triangle':
                img_res = img <= filt.threshold_triangle(img)
            
            img_res = img_as_ubyte(img_res)

            index1 = np.where(img_res==0)
            index2 = np.where(img_res==255)
            if len(index1[0])>len(index2[0]):
                img_res[index1]=255
                img_res[index2]=0
                if invert:
                    img_res[index1]=0
                    img_res[index2]=255
            else:
                img_res[index2]=255
                img_res[index1]=0
                if invert:
                    img_res[index2]=0
                    img_res[index1]=255

            if self.input_settings['FiberData']:
                concomp = measure.label(img_res)
                props = measure.regionprops(concomp)
                for k in range(len(props)):
                    if props[k].eccentricity>0.97 or props[k].area<100:
                        label = props[k].label
                        img_res[concomp==label]=0

            if in_memory:
                self.segmented[i,:,:] = img_res
            else:
                name = str(output_tif_dir)+"/Slice"+"{0:0>4}".format(j)+".tif"
                misc.imsave(name,img_res)
                self.segmented_filenames.append(name)

            print(colored('Finished image '+str(j),'yellow','on_grey', attrs=['bold']))
            j = j+1

