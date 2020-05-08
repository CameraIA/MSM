"""
.. module:: SegmentationAlgorithm
   :synopsis: This module defines a segmentation algorithm

.. moduleauthor:: T. Perciano

"""

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

from abc import ABC, abstractmethod

class SegmentationAlgorithm(ABC):
    """ An abstract class of a segmentation algorithm

    :param image:               Input image
    :type image:                ndimage or list
    :param input_settings:      Dictionary with input settings
    :type input_settings:       dic
    :param preproc_settings:    Dictionaty with preprocessing settings
    :type preproc_settings:     dic
    :param seg_settings:        Dictionaty with segmentation settings
    :type seg_settings:         dic
    :returns:                   Segmented image
    :rtype:                     ndimage

    """

    def __init__(self, image, input_settings, preproc_settings,seg_settings,run_number):
        super().__init__()
        self.image = image
        self.input_settings = input_settings
        self.preproc_settings = preproc_settings
        self.seg_settings = seg_settings
        self.run_number = run_number
        self.segmented = None
        self.segmented_filenames = list()

    def getResult(self):
        return self.segmented
        
    def getResultFilenames(self):
        return self.segmented_filenames

    @abstractmethod
    def segment(self):
        pass
