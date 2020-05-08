"""
    .. module:: Util
    :synopsis: This module is responsible for providing auxiliary functions for MSMcam pipeline

    .. moduleauthor:: T. Perciano

    """

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

import re
import os
import shutil
from scipy import ndimage as ndi
import numpy as np
from colorama import init
from termcolor import colored
import sys
if sys.version_info[0] < 3:
    from ConfigParser import ConfigParser
else:
    from configparser import ConfigParser
    
class MSMConfigParser(ConfigParser, object):
    def __init__(self, *args, **kwargs):
        super(MSMConfigParser, self).__init__(*args, **kwargs)
        
    def getvals(self, section):
        if not section in self.sections():
            raise ValueError('section not in config')
        res = dict()
        for k, v in self.items(section):
            res[k] = self.uni2val(v)
        return res
            
    @staticmethod
    def uni2val(a):
        if a == 'True':
            return True
        elif a == 'False':
            return False
        elif a == 'None':
            return None
        elif re.match('\d+\.\d+', a):
            return float(a)
        elif re.match('\d+', a):
            return int(a)
        else:
            return str(a)

def create_dir(name):
    """ Creates new directory named 'name'

    :param name:    Desired name for directory
    :type name:     str

    """
    try:
        os.makedirs(str(name))
    except OSError:
        shutil.rmtree(str(name))
        os.makedirs(str(name))

def upscale(input_image, scale):
    """ Upscale image

    :param input_image: Input image
    :type input_image:  ndimage
    :param scale:       Scale factor
    :type scale:        int
    :returns:           Scaled image
    :rtype:             ndimage

    """
    out = ndi.zoom(input_image,scale)
    index = out>255
    out[index] = 255
    index = out<0
    out[index] = 0
    #out = out.astype(np.uint8)
    return(out)

def listdir_fullpath(d):
    """ Creates array of filenames with full path given input directory

    :param dir:     Input directory
    :type dir:      str
    :returns:       Array of file paths with full names
    :rtype:         array

    """

    return [os.path.join(d, f) for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

def print_error(msg):
    """ Prints colored error message

    :param msg:     Message
    :type:          str

    """
    print(colored(msg,'red','on_grey',attrs=['bold']))
    os._exit(0)

def check_settings(settings):
    """ Treats error in the parameter file .ini

    :param settings:    Dictionary with settings
    :type settings:     dic

    """
    input_settings = settings.getvals('Input')
    try:
        in_type = input_settings["InputType"]
        if in_type!=0 and in_type!=1:
            print_error('Error in input type, please check your settings file!')
    except:
        print_error('Error in input type, please check your settings file!')

    in_dir = input_settings["InputDir"]
    if os.path.exists(in_dir) is False and in_type==0:
        print_error('Input directory does not exist. If you informed a file name because image type is a single image stack, please change InputType to 1. Otherwise, please check informed InputDir!')
    in_dir_split = in_dir.split("/")
    last = in_dir_split[-1]
    if in_type==1 and "." not in last:
        print_error('If image type is a single image stack, please inform full path including file name in InputDir.')
    if in_type==0 and "." in last:
        print_error('If image type is a sequence of files, please inform only full path in InputDir.')
    try:
        gt_bool = input_settings["GroundTruthExists"]
        gt_dir = input_settings["GroundTruthDir"]
        if gt_bool is not True and gt_bool is not False:
            print_error('Error in GroundTruthExists, please check your settings file!')
        else:
            if gt_bool:
                if os.path.exists(gt_dir) is False and in_type==0:
                    print_error('Groundtruth directory does not exist. If you informed a file name because image type is a single image stack, please change InputType to 1. Otherwise, please check informed GroundTruthDir!')
                gt_dir_split = gt_dir.split("/")
                last = gt_dir_split[-1]
                if in_type==1 and "." not in last:
                    print_error('If image type is a single image stack, please inform full path including file name in GroundTruthDir.')
                if in_type==0 and "." in last:
                    print_error('If image type is a sequence of files, please inform only full path in GroundTruthDir.')
    except:
        print_error('Error in GroundTruthExists, please check your settings file!')

    try:
        fd_bool = input_settings["FiberData"]
        if fd_bool is not True and fd_bool is not False:
            print_error('Error in FiberData, please check your settings file!')
    except:
        print_error('Error in FiberData, please check your settings file!')

    try:
        first_slice = int(input_settings["FirstSlice"])
        if first_slice<0:
            print_error('Error in FirstSlice, please check your settings file!')
    except:
        print_error('Error in FirstSlice, please check your settings file!')
    
    try:
        last_slice = int(input_settings["LastSlice"])
        if last_slice<=0 or last_slice<=first_slice:
            print_error('Error in LastSlice, please check your settings file!')
    except:
        print_error('Error in LastSlice, please check your settings file!')

    if in_type==0:
        number_of_files = len(listdir_fullpath(in_dir))
        if last_slice>number_of_files:
            print_error('Error in LastSlice, last slice number does not match the number of files available, please check your settings file!')

    preproc_settings = settings.getvals("Preprocess")
    try:
        d_scale = preproc_settings["DownsampleScale"]
        if d_scale<=0:
            print_error('Error in DownsampleScale, please check your settings file!')
    except:
        print_error('Error in DownsampleScale, please check your settings file!')

    try:
        use_med = preproc_settings["UseMedian"]
        if use_med is not True and use_med is not False:
            print_error('Error in UseMedian, please check your settings file!')
    except:
        print_error('Error in UseMedian, please check your settings file!')

    try:
        mm_size = preproc_settings["MedianMaskSize"]
        if mm_size<=1:
            print_error('Error in MedianMaskSize, please check your settings file!')
    except:
        print_error('Error in MedianMaskSize, please check your settings file!')

    try:
        use_bil = preproc_settings["UseBilateral"]
        if use_bil is not True and use_bil is not False:
            print_error('Error in UseBilateral, please check your settings file!')
    except:
        print_error('Error in UseBilateral, please check your settings file!')

    try:
        bil_sigcolor = preproc_settings["BilateralSigmaColor"]
        if bil_sigcolor<0 or bil_sigcolor>1:
            print_error('Error in BilateralSigmaColor, please check your settings file!')
    except:
        print_error('Error in BilateralSigmaColor, please check your settings file!')

    try:
        bil_sigspatial = preproc_settings["BilateralSigmaSpatial"]
        if bil_sigspatial<0:
            print_error('Error in BilateralSigmaSpatial, please check your settings file!')
    except:
        print_error('Error in BilateralSigmaSpatial, please check your settings file!')
    seg_settings = settings.getvals("Segmentation")
    try:
        pmrf = seg_settings["RunPMRF"]
        if pmrf is not True and pmrf is not False:
            print_error('Error in RunPMRF, please check your settings file!')
    except:
        print_error('Error in RunPMRF, please check your settings file!')

    try:
        srm = seg_settings["RunSRM"]
        if srm is not True and srm is not False:
            print_error('Error in RunSRM, please check your settings file!')
    except:
        print_error('Error in RunSRM, please check your settings file!')

    try:
        kmeans = seg_settings["RunKMEANS"]
        if kmeans is not True and kmeans is not False:
            print_error('Error in RunKMEANS, please check your settings file!')
    except:
        print_error('Error in RunKMEANS, please check your settings file!')
