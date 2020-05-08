"""
    .. module:: PMRF
    :synopsis: This module is responsible for image segmentation using pmrf algorithm
    
    .. moduleauthor:: T. Perciano
    
    """

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

import ctypes

from ..srm.pysrm.srm import segment_aux
import numpy as np
from scipy import misc
import subprocess
import os
from ...util.util import create_dir, listdir_fullpath, upscale
from skimage import io, img_as_float
from colorama import init
from termcolor import colored
from ..SegmentationAlgorithm import SegmentationAlgorithm

class PMRF(SegmentationAlgorithm):
    def __init__(self, image, input_settings, preproc_settings,seg_settings,run_number):
        super().__init__(image, input_settings, preproc_settings,seg_settings,run_number)

    def segment(self):
        """ Main segmentation process using PMRF algorithm
        
        """
        input_type = self.input_settings["InputType"]
        input_dir = self.input_settings["InputDir"]
        if input_type==1:
            input_dir_split = input_dir.split("/")
            input_dir = "/".join(input_dir_split[:-1])
            input_file = input_dir
        output_dat_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number),"res/pmrf/dat")
        output_tif_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number), "res/pmrf/tif")
        input_dat_dir = os.path.join(input_dir, "msmcam_run_"+str(self.run_number),"dat")
        overseg_tif_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number), "preproc/overseg/tif")
        overseg_dat_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number), "preproc/overseg/dat")
        in_memory = self.input_settings["InMemory"]
        first_slice = self.input_settings["FirstSlice"]
        last_slice = self.input_settings["LastSlice"]
        num_slices = last_slice-first_slice
        if in_memory:
            self.segmented = np.zeros(self.image.shape,dtype=np.uint8)
        
        create_dir(overseg_tif_dir)
        create_dir(overseg_dat_dir)
        
        # Prepare oversegmentations using srm
        #if in_memory:
        #    self.image = upscale(self.image, self.preproc_settings["DownsampleScale"])
        j = first_slice
        for i in range(num_slices):
            if in_memory:
                img_slice = self.image[i,:,:]*1.0
            else:
                img_slice = io.imread(str(self.image[i]))*1.0
            avg_out, lbl_out = segment_aux(img_slice, q=40)
            name = str(overseg_tif_dir)+"/Slice"+"{0:0>4}".format(j)+".tif"
            misc.imsave(name,avg_out)
            name = str(overseg_dat_dir)+"/Slice"+"{0:0>4}".format(j)+".dat"
            avg_out = avg_out.astype(np.uint8)
            with open(name, 'wb') as f:
                avg_out.tofile(f)
            j = j+1

        ###########################
        
        filenames = sorted(listdir_fullpath(str(input_dat_dir)))
        
        create_dir(output_dat_dir)
        create_dir(output_tif_dir)
       
        multiphase = self.seg_settings["Multiphase"]
        invert = self.seg_settings["Invert"]
        
        if multiphase:
            num_labels = self.seg_settings["NumClustersPMRF"]
        else:
            num_labels = 2
        
        if in_memory:
            z, y, x = self.image.shape
        else:
            y, x = io.imread(str(self.image[0])).shape
            z = num_slices
        j = first_slice
        for i in range(len(filenames)):
            input_name = filenames[i]
            output_dat_name = str(output_dat_dir)+"/Slice"+"{0:0>4}".format(j)+".dat"
            output_tif_name = str(output_tif_dir)+"/Slice"+"{0:0>4}".format(j)+".tif"
            overseg_name = str(overseg_dat_dir) + "/Slice"+"{0:0>4}".format(j)+".dat"       
                    
            print('python '+os.path.dirname(os.path.abspath(__file__))+'/MPI_PMRF.py -i '+input_name+ ' -s '+ overseg_name+ ' -o '+output_dat_name+' -b '+ str(x)+','+str(y)+','+str(z)+ ' -e 10 -m 10 -l '+str(num_labels))
            
            p = subprocess.Popen(['python '+os.path.dirname(os.path.abspath(__file__))+'/MPI_PMRF.py -i '+input_name+ ' -s '+ overseg_name+ ' -o '+output_dat_name+' -b '+ str(x)+','+str(y)+','+str(z)+ ' -e 10 -m 10 -l '+str(num_labels)], shell=True)
            (output, err) = p.communicate()
            p_status = p.wait()
                    
            #Saving .tif from .dat
            with open(output_dat_name, 'rb') as f:
                img_res = np.fromfile(f, dtype=np.uint8)
                img_res.shape = (y,x)
                if num_labels==2 or multiphase is False:
                    index1 = np.where(img_res==30)
                    index2 = np.where(img_res==60)
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
                misc.imsave(output_tif_name,img_res)
                if in_memory:
                    self.segmented[i,:,:] = img_res
                
            #print('Finished image '+str(i))
            print(colored('Finished image '+str(j),'yellow','on_grey', attrs=['bold']))
            j =  j+1

        self.segmented_filenames = sorted(listdir_fullpath(str(output_tif_dir)))

