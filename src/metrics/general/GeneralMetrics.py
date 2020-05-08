"""
.. module:: GeneralMetrics
   :synopsis: This module is responsible for calculating the general metrics. Metrics are calculated for the entire 3D stack.

.. moduleauthor:: T. Perciano

"""

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

import numpy as np
import scipy as sp
from ...util.util import create_dir
from sklearn.metrics import confusion_matrix
from ...inout.ImageReader import ImageReader
from skimage import measure, io
from colorama import init
from termcolor import colored

class GeneralMetrics:
    """ A class containing general metrics for segmented images

    :param seg_image:		Segmented image
    :type seg_image:		ndimage
    :param input_settings:	Dictionary with input settings

    """

    def __init__(self, segmented, input_settings,run_number):
        """ Initialize GeneralMetrics

        """
        #### self.seg_image = seg_image
        self.segmented = segmented
        """ Segmented image """
        self.gt_image = None
        self.gt_filenames = None
        """ Groundtruth image """
        self.conf_mat_values = None
        """ Confusion matrix values """
        self.input_type = input_settings["InputType"]
        """ Image input type (sequence or stack) """
        self.input_dir = input_settings["InputDir"]
        """ Name of the input directory """
        self.gt_dir = input_settings["GroundTruthDir"]
        """ Name of the groundtruth directory """
        self.first_slice = int(input_settings["FirstSlice"])
        """ First slice """
        self.last_slice = int(input_settings["LastSlice"])
        """ Last slice """
        self.num_slices = self.last_slice-self.first_slice
        """ Number of slices """
        self.in_memory = input_settings["InMemory"]
        self.run_number = run_number

        self.total_ones_gt = 0
        self.total_zeros_gt = 0
        self.total_ones_seg = 0
        self.total_zeros_seg = 0

    def getMetrics(self):
        """ Returns the metrics

        :returns:	List with metrics values
        :rtype:     list

        """

        return(self.metrics)

    def getGT_filenames(self):
        """ Returns groundtruth filenames """
        return(self.gt_filenames)

    def getGT_image(self):
        """ Returns groundtruth """
        return(self.gt_image)

    def read_groundtruth(self):
        """ Reads groundtruth data """
        reader = ImageReader(self.gt_dir, self.input_type,self.first_slice, self.last_slice,False,self.in_memory,self.run_number)
        reader.read()
        if self.in_memory:
        	self.gt_image = reader.getImage()
        else:
        	self.gt_filenames = reader.getImageFilenames()

    def calc_conf_matrix(self):
        """ Calculates confusion matrix using groundtruth """
        tn_global = 0
        fp_global = 0
        fn_global = 0
        tp_global = 0
        if self.in_memory:
            z, y, x = self.gt_image.shape
        else:
            y, x = io.imread(str(self.gt_filenames[0])).shape
            z = self.num_slices
        j = self.first_slice
        for i in range(z):

            if self.in_memory:
                gt_image = self.gt_image[i,:,:]
                seg_image = self.segmented[i,:,:]
            else:
                gt_image = io.imread(str(self.gt_filenames[i]))
                seg_image = io.imread(str(self.segmented[i]))

            index1 = np.where(gt_image==0)
            index2 = np.where(gt_image==255)

            if len(index1[0])>len(index2[0]):
                gt_image[index1]=1
                gt_image[index2]=0
            else:
                gt_image[index2]=1
                gt_image[index1]=0

            index1 = np.where(seg_image==0)
            index2 = np.where(seg_image==255)
            if len(index1[0])>len(index2[0]):
                seg_image[index1]=1
                seg_image[index2]=0
            else:
                seg_image[index2]=1
                seg_image[index1]=0

            gt_array = gt_image.flatten()
            seg_array = seg_image.flatten()

            self.total_zeros_gt = self.total_zeros_gt + float(sp.sum(gt_array==0))
            self.total_ones_gt = self.total_ones_gt + float(sp.sum(gt_array==1))
            self.total_zeros_seg = self.total_zeros_seg + float(sp.sum(seg_array==0))
            self.total_ones_seg = self.total_ones_seg + float(sp.sum(seg_array==1))

            conf_mat = confusion_matrix(gt_array, seg_array)
            tn, fp, fn, tp = conf_mat.ravel()
            tn_global = tn_global + tn
            fp_global = fp_global + fp
            fn_global = fn_global + fn
            tp_global = tp_global + tp

            print(colored('Finished image '+str(j),'yellow','on_grey', attrs=['bold']))
            j = j+1

        self.conf_mat_values = (tn_global,fp_global,fn_global,tp_global)

    def calc_porosity(self):
        """ Calculates porosity of both segmented data and groundtruth

        :returns:   List with porosity values for segmented and groundtruth data
        :rtype:     list

        """
        #### Vp = float(sp.sum(self.seg_array == 0))
        Vp = self.total_zeros_seg
        #### Vs = float(sp.sum(self.seg_array == 1))
        Vs = self.total_ones_seg
        seg_porosity = round(Vp/(Vs + Vp),3)
        #### Vp = float(sp.sum(self.gt_array == 0))
        Vp = self.total_zeros_gt
        #### Vs = float(sp.sum(self.gt_array == 1))
        Vs = self.total_ones_gt
        gt_porosity = round(Vp/(Vs + Vp),3)
        return(seg_porosity,gt_porosity)

    def calc_volume(self):
        """ Calculates volume of both segmented data and groundtruth

        :returns:   List with volume values for segmented and groundtruth data
        :rtype:     list

        """
        #### seg_volume = float(sp.sum(self.seg_array == 1))
        seg_volume = self.total_ones_seg
        #### gt_volume = float(sp.sum(self.gt_array == 1))
        gt_volume = self.total_ones_gt
        return(seg_volume,gt_volume)

    def process(self):
        """ Full pipeline to calculate general metrics """
        self.read_groundtruth()
        self.calc_conf_matrix()
        tn,fp,fn,tp = self.conf_mat_values
        tn = float(tn)
        fp = float(fp)
        fn = float(fn)
        tp = float(tp)
        precision = round(tp/(tp+fp),3)
        recall = round(tp/(tp+fn),3)
        accuracy = round((tp+tn)/(tp+tn+fp+fn),3)
        specificity = round(tp/(tn+fp),3)
        F = round(2*((precision*recall)/(precision+recall)),3)
        seg_porosity, gt_porosity = self.calc_porosity()
        seg_volume, gt_volume = self.calc_volume()
        self.metrics = (precision,recall,accuracy,specificity,F,seg_porosity,gt_porosity,seg_volume,gt_volume)
