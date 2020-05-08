"""
.. module:: MetricsWriter
   :synopsis: This module is responsible for saving general metrics to a spreadsheet

.. moduleauthor:: T. Perciano

"""

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"


import os
from ...util.util import create_dir
import pandas as pd
from ...vis.Visualizer import Visualizer

class MetricsWriter:
    """ A class to create a table with general metrics for different algorithms

    :param algorithms:		List of segmentation algorithms
    :type algorithms:       list
    :param metrics:         List of metrics
    :type metrics:          list
    :param input_dir:       Name of input directory
    :type input_dir:        str
    :param vis_settings:    Dictionary with visualization settings
    :type vis_settings:     dic

    """

    def __init__(self, algorithms, metrics, input_dir, input_type, vis_settings, in_memory, run_number):
        """ Initialize MetricsWriter

        """
        self.input_dir_back = input_dir
        self.input_type = input_type
        if self.input_type==1:
            input_dir_split = input_dir.split("/")
            input_dir = "/".join(input_dir_split[:-1])
            input_file = input_dir
        self.algorithms = algorithms
        """ List of segmentation algorithms """
        self.metrics = metrics
        """ List of metric values """
        self.run_number = run_number
        self.output_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number),'metrics/general')
        """ Output directory """
        self.input_dir = input_dir
        """ Input directory """
        create_dir(self.output_dir)
        self.output_file = os.path.join(self.output_dir,'metrics.xlsx')
        """ File name for spreadsheet output """
        self.vis_settings = vis_settings
        """ Dictionary with visualization settings """
        self.in_memory = in_memory

    def createTable(self):
        """ Function to create table with general metrics and visualize graphics """
        df = pd.DataFrame(data=self.metrics,index=self.algorithms,columns=['Precision','Recall','Accuracy','Specificity','F','SegPorosity','GtPorosity','SegVolume','GtVolume'])
        writer = pd.ExcelWriter(self.output_file, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='Sheet1')
        writer.save()

        # Vis
        vis = Visualizer(self.input_dir_back, self.input_type, self.vis_settings, self.in_memory, self.run_number)
        vis.viewGeneralMetrics(df, self.output_dir)
