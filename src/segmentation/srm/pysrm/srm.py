"""
    .. module:: SRM
    :synopsis: This module is responsible for image segmentation using srm algorithm
    
    .. moduleauthor:: D. Ushizima
    
    """

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "D. Ushizima"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

import os
import numpy as np
from ctypes import cdll, c_float, c_int32, c_ubyte, c_bool, POINTER
from skimage import img_as_float, img_as_ubyte, segmentation, filters
from skimage import morphology, measure
from scipy import misc
import platform
import matplotlib.pyplot as plt
from skimage import io
from colorama import init
from termcolor import colored
from ....util.util import create_dir, upscale, listdir_fullpath
from ...SegmentationAlgorithm import SegmentationAlgorithm

if platform.system() == "Darwin":
    srm_lib = cdll.LoadLibrary(os.path.dirname(os.path.abspath(__file__))+'/libsrm.dylib')
else:
    srm_lib = cdll.LoadLibrary(os.path.dirname(os.path.abspath(__file__))+'/libsrm.so')
srm_fn = srm_lib.srm_c

def segment_aux(img, q=25, get_labels=True, get_average=True):
    """
    Performs image segmentation using the SRM algorithm 
    
    :param img:             input image
    :type img:              ndarray
    :param q:               complexity of the assumed distributions (affects the number of unique colors)
    :type q:                int
    :param get_labels:      compute and return the labeled image
    :type get_labels:       bool
    :param get_average:     compute and return the image with region's average colors
    :type get_average:      bool
    :returns:               tuple of output images
    :rtype:                 (ndarray, ndarray)
    """
    if not img.dtype == np.ubyte:
        min = img.min()
        max = img.max()
        img = (img - min) / (max - min) * 255
        img = img.astype(np.ubyte)
    c_q = c_float(q)
    c_img = img.ctypes.data_as(POINTER(c_ubyte))
    c_w = c_int32(img.shape[1])
    c_h = c_int32(img.shape[0])
    avg_out = np.zeros(img.shape, dtype=np.float32)
    lbl_out = np.zeros(img.shape, dtype=np.int32)
    c_avg_out = avg_out.ctypes.data_as(POINTER(c_float))
    c_lbl_out = lbl_out.ctypes.data_as(POINTER(c_int32))
    srm_fn(c_q, c_img, c_w, c_h, c_bool(get_average), c_bool(get_labels), c_avg_out, c_lbl_out)
    return avg_out, lbl_out

class SRM(SegmentationAlgorithm):

    def __init__(self, image, input_settings, preproc_settings,seg_settings,run_number):
        super().__init__(image, input_settings, preproc_settings,seg_settings,run_number)

        self.thresh_alg = None

    def onclick(self,event):
        '''
        Event handler for button_press_event
        @param event MouseEvent
        '''
        self.thresh_alg = event.inaxes.get_title()
        plt.close()
        
    def segment(self):
        """ Main segmentation processing using srm algorithm 
        
        """      
        input_type = self.input_settings["InputType"]
        input_dir = self.input_settings["InputDir"]
        if input_type==1:
            input_dir_split = input_dir.split("/")
            input_dir = "/".join(input_dir_split[:-1])
            input_file = input_dir
        
        masked = self.input_settings["Masked"]
        invert = self.seg_settings["Invert"]
        qsrm = self.seg_settings["QSRM"]
        in_memory = self.input_settings["InMemory"]
        first_slice = self.input_settings["FirstSlice"]
        last_slice = self.input_settings["LastSlice"]
        num_slices = last_slice-first_slice

        if in_memory:
            self.segmented = np.zeros(self.image.shape,dtype=np.uint8)

        if not in_memory:
            output_tif_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number), "res/srm/tif")
            create_dir(output_tif_dir)

        minSize = 50
        multiphase = self.seg_settings["Multiphase"]
        
        if multiphase is False:
	        if masked:
	            n_colors = 3
	        else:
	            n_colors = 2
          
        #if in_memory:
        #    self.image = upscale(self.image, self.preproc_settings["DownsampleScale"])

        j = first_slice
        for i in range(num_slices):
            if in_memory:
                img2 = img_as_float(self.image[i,:,:])
            else:
                img2 = img_as_float(io.imread(str(self.image[i])))

            avg_out, lbl_out = segment_aux(img2, q=qsrm)
            if multiphase is False:
                if masked is False:
                    mask = avg_out > filters.threshold_otsu(avg_out)
                    img_res = img_as_ubyte(mask)
                    index1 = np.where(img_res==0)
                    index2 = np.where(img_res==255)
                else:
                    if i==0:
                        fig, ax = filters.try_all_threshold(avg_out, figsize=(10, 6), verbose=False)
                        cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
                        plt.suptitle('Please click on the best result figure', fontsize=15)
                        plt.show()
                    if self.thresh_alg=="Li":
                        mask = avg_out > filters.threshold_li(avg_out)
                    if self.thresh_alg=="Triangle":
                        mask = avg_out > filters.threshold_triangle(avg_out)
                    if self.thresh_alg=="Mean":
                        mask = avg_out > filters.threshold_mean(avg_out)
                    if self.thresh_alg=="Isodata":
                        mask = avg_out > filters.threshold_isodata(avg_out)
                    if self.thresh_alg=="Otsu":
                        mask = avg_out > filters.threshold_otsu(avg_out)
                    if self.thresh_alg=="Yen":
                        mask = avg_out > filters.threshold_yen(avg_out)
                    if self.thresh_alg=="Minimum":
                        mask = avg_out > filters.threshold_minimum(avg_out)
                        
                    img_res = img_as_ubyte(mask)
                    label_background = img_res[0,0]
                    print(label_background)
                    count = 0
                    for label in range(n_colors):
                        if label != label_background:
                            if count==0:
                                index1 = np.where(img_res==label)
                            if count==1:
                                index2 = np.where(img_res==label)
                            count = count+1
                    img_res[img_res==label_background]=0    
                    #index1 = np.where(img_res==0)
                    #index2 = np.where(img_res==255)
                               
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
            else:
                if in_memory:
                    self.segmented[i,:,:] = avg_out
                else:
                    name = str(output_tif_dir)+"/Slice"+"{0:0>4}".format(j)+".tif"
                    misc.imsave(name,avg_out)
                    self.segmented_filenames.append(name)
            print(colored('Finished image '+str(j),'yellow','on_grey', attrs=['bold']))
            j = j+1
            
