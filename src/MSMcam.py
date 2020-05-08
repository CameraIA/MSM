#!/usr/bin/env python
"""MSMcam main file

This is the CAMERA Materials Segmentation & Metrics (MSMcam) package as described in the Journal of Synchrotron Radiation article: "Insight into 3D micro-CT data: exploring segmentation algorithms through performance metrics" by Perciano, Ushizima, Krishnan, Parkinson, Larson, Pelt, Bethel, Zok and Sethian, although the codes also count with contributions from J. Mike Macneil.

[Full paper] http://onlinelibrary.wiley.com/doi/10.1107/S1600577517010955/abstract

"""

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "CAMERA Image Processing"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

import argparse
import sys
import numpy as np
import configparser
import os

from .inout.ImageReader import ImageReader
from .preprocessing.PreProcessor import PreProcessor
from .postprocessing.PostProcessor import PostProcessor
from .util.util import create_dir, upscale, listdir_fullpath, check_settings
from .util.util import MSMConfigParser

from .segmentation.srm.pysrm.srm import SRM
from .segmentation.pmrf.pmrf import PMRF
from .segmentation.kmeans.kmeans import KMEANS
from .segmentation.threshold.threshold import Threshold
from .metrics.general.GeneralMetrics import GeneralMetrics
from .metrics.general.MetricsWriter import MetricsWriter
from .metrics.fibers.FiberMetrics import FiberMetrics, make_plots, write_pdf_report
from .vis.Visualizer import Visualizer

from colorama import init
from termcolor import colored
from glob import glob
from shutil import copyfile

def parse_cmdl_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='CAM_MSM')
    parser.add_argument("-c", "--config", type=str, default=None, help="Configuration file")
    return parser.parse_args()

def load_settings(file_name):
    """
    Helper function to parse the INI file.

    :param file_name:   Settings file name
    :type file_name:    str
    :returns:           Dictionary with settings from .ini file
    :rtype:             dic
    """
    settings = MSMConfigParser()
    settings.optionxform = str
    settings.read(file_name)
    return(settings)

def main():
    """MSMcam main function"""
    import time
    t_total = time.time()
    image_collection_filenames = { "Orig" : None, "GT" : None, "Filtered" : None, "SegPMRF" : None, "SegSRM" : None, "SegKMEANS" : None , "SegThresh" : None}
    image_collection = { "Orig" : None, "GT" : None, "Filtered" : None, "SegPMRF" : None, "SegSRM" : None, "SegKMEANS" : None , "SegThresh" : None}
    init()
    print(colored('Starting MSMcam...','green','on_grey', attrs=['bold']))
    args = parse_cmdl_args()
    settings = load_settings(args.config)
    check_settings(settings)

    input_settings = settings.getvals('Input')

    input_type = input_settings["InputType"]
    input_dir = input_settings["InputDir"]
    if input_type==1:
        input_dir_split = input_dir.split("/")
        input_dir = "/".join(input_dir_split[:-1])+"/"
    paths_list = glob(input_dir+"*/")
    temp = list()
    for i in range(len(paths_list)):
        if "msmcam_run" in paths_list[i]:
            temp.append(paths_list[i])
    paths_list = sorted(temp)
    if len(paths_list)>0:
        last_run = paths_list[-1]
        last_run = last_run.split("/")[-2]
        run_number = int(last_run.split("_")[-1]) + 1
    else:
        run_number = 0
    create_dir(input_dir+"/msmcam_run_"+str(run_number))
    copyfile(args.config,input_dir+"/msmcam_run_"+str(run_number)+"/parameters.ini")

    first_slice = input_settings["FirstSlice"]
    last_slice = input_settings["LastSlice"]
    num_slices = last_slice-first_slice
    has_gt = input_settings["GroundTruthExists"]
    if has_gt:
        gt_dir = input_settings["GroundTruthDir"]
    fiber_data = input_settings["FiberData"]
    preproc_settings = settings.getvals('Preprocess')
    seg_settings = settings.getvals('Segmentation')
    postproc_settings = settings.getvals('Postprocess')
    vis_settings = settings.getvals("Visualization")
    in_memory = input_settings["InMemory"]

    print(colored('Reading images...','green','on_grey', attrs=['bold']))

    reader = ImageReader(input_settings["InputDir"],int(input_settings["InputType"]),first_slice,last_slice,True,in_memory,run_number)
    t_read = time.time()
    reader.read()
    t_read = time.time() - t_read
    if in_memory:
        input_image = reader.getImage()
        image_collection["Orig"] = input_image
    input_filenames = reader.getImageFilenames()
    image_collection_filenames["Orig"] = input_filenames
    print(colored('Finished reading images!','green','on_grey', attrs=['bold']))

    # Preprocessing
    print(colored('Preprocessing...','green','on_grey', attrs=['bold']))
    if in_memory:
        preprocessor = PreProcessor(input_image, input_settings,preproc_settings,run_number)
    else:
        preprocessor = PreProcessor(input_filenames, input_settings,preproc_settings,run_number)
    t_filtering = time.time()
    preprocessor.process()
    t_filtering = time.time() - t_filtering
    input_filtered, filtered_filenames = preprocessor.getFiltered()
    image_collection_filenames["Filtered"] = filtered_filenames
    if in_memory:
        image_collection["Filtered"] = input_filtered
    print(colored('Finished preprocessing!','green','on_grey', attrs=['bold']))

    print(colored('Running segmentation algorithms...','green','on_grey', attrs=['bold']))
    # Run PMRF
    if seg_settings["RunPMRF"]:
        print(colored('PMRF','yellow','on_grey', attrs=['bold']))
        if in_memory:
            pmrf_alg = PMRF(input_filtered, input_settings, preproc_settings, seg_settings,run_number)
        else:
            pmrf_alg = PMRF(filtered_filenames, input_settings, preproc_settings, seg_settings,run_number)
        t_pmrf = time.time()
        pmrf_alg.segment()
        t_pmrf = time.time() - t_pmrf
        
        res = pmrf_alg.getResult()
        res_files = pmrf_alg.getResultFilenames()
        # Postprocessing
        if postproc_settings["RunPostProc"]:
            if in_memory:
                postprocessor = PostProcessor(res, input_settings, postproc_settings,run_number)
            else:
                postprocessor = PostProcessor(res_files, input_settings, postproc_settings,run_number)
            postprocessor.process()
            res = postprocessor.getResult()
            res_files = postprocessor.getResultFilenames()
        ################
        
        image_collection["SegPMRF"] = res
        image_collection_filenames["SegPMRF"] = res_files

    # Run SRM
    if seg_settings["RunSRM"]:
        print(colored('SRM','yellow','on_grey', attrs=['bold']))
        if in_memory:
            srm_alg = SRM(input_filtered, input_settings, preproc_settings, seg_settings,run_number)
        else:
            srm_alg = SRM(filtered_filenames, input_settings, preproc_settings, seg_settings,run_number)
        t_srm = time.time()
        srm_alg.segment()
        t_srm = time.time() - t_srm
        
        res = srm_alg.getResult()
        res_files = srm_alg.getResultFilenames()
        # Postprocessing
        if postproc_settings["RunPostProc"]:
            if in_memory:
                postprocessor = PostProcessor(res, input_settings, postproc_settings,run_number)
            else:
                postprocessor = PostProcessor(res_files, input_settings, postproc_settings,run_number)
            postprocessor.process()
            res = postprocessor.getResult()
            res_files = postprocessor.getResultFilenames()
        ################
        
        image_collection["SegSRM"] = res
        image_collection_filenames["SegSRM"] = res_files

    # Run KMEANS
    if seg_settings["RunKMEANS"]:
        print(colored('KMEANS','yellow','on_grey', attrs=['bold']))
        if in_memory:
            kmeans_alg = KMEANS(input_filtered, input_settings, preproc_settings, seg_settings,run_number)
        else:
            kmeans_alg = KMEANS(filtered_filenames, input_settings, preproc_settings, seg_settings,run_number)
        t_kmeans = time.time()
        kmeans_alg.segment()
        t_kmeans = time.time() - t_kmeans
        
        res = kmeans_alg.getResult()
        res_files = kmeans_alg.getResultFilenames()
        # Postprocessing
        if postproc_settings["RunPostProc"]:
            if in_memory:
                postprocessor = PostProcessor(res, input_settings, postproc_settings,run_number)
            else:
                postprocessor = PostProcessor(res_files, input_settings, postproc_settings,run_number)
            postprocessor.process()
            res = postprocessor.getResult()
            res_files = postprocessor.getResultFilenames()
        ################
        
        image_collection["SegKMEANS"] = res
        image_collection_filenames["SegKMEANS"] = res_files

    # Run Threshold
    if seg_settings["RunThresh"]:
        print(colored('Threshold','yellow','on_grey', attrs=['bold']))
        if in_memory:
            thresh_alg = Threshold(input_filtered, input_settings, preproc_settings, seg_settings,run_number)
        else:
            thresh_alg = Threshold(filtered_filenames, input_settings, preproc_settings, seg_settings,run_number)
        t_thresh = time.time()
        thresh_alg.segment()
        t_thresh = time.time() - t_thresh
        
        res = thresh_alg.getResult()
        res_files = thresh_alg.getResultFilenames()
        # Postprocessing
        if postproc_settings["RunPostProc"]:
            if in_memory:
                postprocessor = PostProcessor(res, input_settings, postproc_settings,run_number)
            else:
                postprocessor = PostProcessor(res_files, input_settings, postproc_settings,run_number)
            postprocessor.process()
            res = postprocessor.getResult()
            res_files = postprocessor.getResultFilenames()
        ################
        
        image_collection["SegThresh"] = res
        image_collection_filenames["SegThresh"] = res_files

    print(colored('Finished segmenting images!','green','on_grey', attrs=['bold']))

    t_general = time.time()
    if has_gt and seg_settings["Multiphase"] is False:
        # Calculate metrics
        print(colored('Calculating general metrics...','green','on_grey', attrs=['bold']))
        algs = list()
        metrics = list()
        # General
        if image_collection_filenames["SegPMRF"] is not None or image_collection["SegPMRF"] is not None:
            print(colored('PMRF','yellow','on_grey', attrs=['bold']))
            if in_memory:
                gen_metrics_pmrf = GeneralMetrics(image_collection["SegPMRF"],input_settings,run_number)
            else:
                gen_metrics_pmrf = GeneralMetrics(image_collection_filenames["SegPMRF"],input_settings,run_number)
            gen_metrics_pmrf.process()
            pmrf_res = gen_metrics_pmrf.getMetrics()
            algs.append('pmrf')
            metrics.append(pmrf_res)
            image_collection_filenames["GT"] = gen_metrics_pmrf.getGT_filenames()
            if in_memory:
                image_collection["GT"] = gen_metrics_pmrf.getGT_image()
        if image_collection_filenames["SegSRM"] is not None or image_collection["SegSRM"] is not None:
            print(colored('SRM','yellow','on_grey', attrs=['bold']))
            if in_memory:
                gen_metrics_srm = GeneralMetrics(image_collection["SegSRM"],input_settings,run_number)
            else:
                gen_metrics_srm = GeneralMetrics(image_collection_filenames["SegSRM"],input_settings,run_number)
            gen_metrics_srm.process()
            srm_res = gen_metrics_srm.getMetrics()
            algs.append('srm')
            metrics.append(srm_res)
            if image_collection_filenames["GT"] is None:
                image_collection_filenames["GT"] = gen_metrics_srm.getGT_filenames()
            if in_memory:
                image_collection["GT"] = gen_metrics_srm.getGT_image()
        if image_collection_filenames["SegKMEANS"] is not None or image_collection["SegKMEANS"] is not None:
            print(colored('KMEANS','yellow','on_grey', attrs=['bold']))
            if in_memory:
                gen_metrics_kmeans = GeneralMetrics(image_collection["SegKMEANS"],input_settings,run_number)
            else:
                gen_metrics_kmeans = GeneralMetrics(image_collection_filenames["SegKMEANS"],input_settings,run_number)
            gen_metrics_kmeans.process()
            kmeans_res = gen_metrics_kmeans.getMetrics()
            algs.append('kmeans')
            metrics.append(kmeans_res)
            if image_collection_filenames["GT"] is None:
                image_collection_filenames["GT"] = gen_metrics_kmeans.getGT_filenames()
            if in_memory:
                image_collection["GT"] = gen_metrics_kmeans.getGT_image()

        if image_collection_filenames["SegThresh"] is not None or image_collection["SegThresh"] is not None:
            print(colored('Threshold','yellow','on_grey', attrs=['bold']))
            if in_memory:
                gen_metrics_thresh = GeneralMetrics(image_collection["SegThresh"],input_settings,run_number)
            else:
                gen_metrics_thresh = GeneralMetrics(image_collection_filenames["SegThresh"],input_settings,run_number)
            gen_metrics_thresh.process()
            thresh_res = gen_metrics_thresh.getMetrics()
            algs.append('thresh')
            metrics.append(thresh_res)
            if image_collection_filenames["GT"] is None:
                image_collection_filenames["GT"] = gen_metrics_thresh.getGT_filenames()
            if in_memory:
                image_collection["GT"] = gen_metrics_thresh.getGT_image()

        mwriter = MetricsWriter(algs,metrics, input_settings["InputDir"], input_settings["InputType"], vis_settings, in_memory,run_number)
        mwriter.createTable()
        print(colored('Finished general metrics!','green','on_grey', attrs=['bold']))

    else:
        print(colored('No groundtruth, skipping general metrics...','green','on_grey', attrs=['bold']))


    t_general = time.time() - t_general
    
    t_fibers = time.time()
    # Fibers
    if fiber_data and seg_settings["Multiphase"] is False:
        # Calculates metrics specific for fibers
        print(colored('Calculating fiber metrics...','green','on_grey', attrs=['bold']))
        if image_collection_filenames["SegPMRF"] is not None or image_collection["SegPMRF"] is not None:
            if in_memory:
                fib_metrics_pmrf = FiberMetrics(image_collection["SegPMRF"],input_settings,'pmrf',run_number)
            else:
                fib_metrics_pmrf = FiberMetrics(image_collection_filenames["SegPMRF"],input_settings,'pmrf',run_number)
            status = fib_metrics_pmrf.process()
            if status is False:
                print(colored('PMRF: Not able to calculate fiber metrics, please check segmentation result. This might be caused by a masked image (with a total black or white background region. Please try with a cropped region without the background.','red','on_grey',attrs=['bold']))
        if image_collection_filenames["SegSRM"] is not None or image_collection_filenames["SegSRM"] is not None:
            if in_memory:
                fib_metrics_srm = FiberMetrics(image_collection_filenames["SegSRM"],input_settings,'srm',run_number)
            else:
                fib_metrics_srm = FiberMetrics(image_collection_filenames["SegSRM"],input_settings,'srm',run_number)
            status = fib_metrics_srm.process()
            if status is False:
                print(colored('SRM: Not able to calculate fiber metrics, please check segmentation result. This might be caused by a masked image (with a total black or white background region. Please try with a cropped region without the background.','red','on_grey',attrs=['bold']))
        if image_collection_filenames["SegKMEANS"] is not None or image_collection_filenames["SegKMEANS"] is not None:
            if in_memory:
                fib_metrics_kmeans = FiberMetrics(image_collection_filenames["SegKMEANS"],input_settings,'kmeans',run_number)
            else:
                fib_metrics_kmeans = FiberMetrics(image_collection_filenames["SegKMEANS"],input_settings,'kmeans',run_number)
            status = fib_metrics_kmeans.process()
            if status is False:
                print(colored('KMEANS: Not able to calculate fiber metrics, please check segmentation result. This might be caused by a masked image (with a total black or white background region. Please try with a cropped region without the background.','red','on_grey',attrs=['bold']))
        if image_collection_filenames["SegThresh"] is not None or image_collection["SegThresh"] is not None:
            if in_memory:
                fib_metrics_thresh = FiberMetrics(image_collection["SegThresh"],input_settings,'thresh',run_number)
            else:
                fib_metrics_thresh = FiberMetrics(image_collection_filenames["SegThresh"],input_settings,'thresh',run_number)
            status = fib_metrics_thresh.process()
            if status is False:
                print(colored('Threshold: Not able to calculate fiber metrics, please check segmentation result. This might be caused by a masked image (with a total black or white background region. Please try with a cropped region without the background.','red','on_grey',attrs=['bold']))
        if has_gt:
            if in_memory:
                fib_metrics_gt = FiberMetrics(image_collection["GT"],input_settings,'gt',run_number)
            else:
                fib_metrics_gt = FiberMetrics(image_collection_filenames["GT"],input_settings,'gt',run_number)
            status = fib_metrics_gt.process()

        print(colored('Ploting graphics...','green','on_grey', attrs=['bold']))
        make_plots(input_settings,run_number)
        print(colored('Writing dashboard...','green','on_grey', attrs=['bold']))
        write_pdf_report(input_settings,run_number)
        print(colored('Finished fibers metrics!','green','on_grey', attrs=['bold']))

    else:
        print(colored('No fiber data, skipping fiber metrics...','green','on_grey', attrs=['bold']))

    t_fibers = time.time()- t_fibers
    # Visualize segmentation results
    print(colored('Starting visualization routine...','green','on_grey', attrs=['bold']))
    if seg_settings["Multiphase"] is False:
        vis = Visualizer(input_settings["InputDir"],input_settings["InputType"], vis_settings, input_settings["InMemory"],run_number)
        overlay, overlay_filenames = vis.viewSegmentations(image_collection, image_collection_filenames)
        if in_memory is False:
            overlay_gif = vis.create_gif()

    # Printing all the output locations
    input_type = input_settings["InputType"]
    input_dir = input_settings["InputDir"]
    if input_type==1:
        input_dir_split = input_dir.split("/")
        input_dir = "/".join(input_dir_split[:-1])+"/"
        
    t_total = time.time() - t_total
    
    print(colored("***********************************************************************************************************************************",'blue','on_white',attrs=['bold']))
    print(colored("All results from the MSMcam pipeline are the following:",'blue','on_white',attrs=['bold']))
    print(colored("Input data is in: "+ input_dir,'blue','on_white',attrs=['bold']))
    print(colored("Segmentation results are in: "+input_dir+"msmcam_run_"+str(run_number)+"/res/", 'blue','on_white',attrs=['bold']))
    print(colored("Results of general and fiber metrics are in: "+input_dir+"msmcam_run_"+str(run_number)+"/metrics/",'blue','on_white',attrs=['bold']))
    print(colored("Results of the preprocessing step are in: "+input_dir+"msmcam_run_"+str(run_number)+"/preproc/",'blue','on_white',attrs=['bold']))
    print(colored("Results of the postprocessing step are in: "+input_dir+"msmcam_run_"+str(run_number)+"/postproc/",'blue','on_white',attrs=['bold']))
    print(colored("Visualization files are in: "+input_dir+"msmcam_run_"+str(run_number)+"/vis/",'blue','on_white',attrs=['bold']))
    print(colored("Elapsed time to read images: "+"{:.3f}".format(t_read)+" seconds.",'blue','on_white',attrs=['bold']))
    print(colored("Elapsed time to filter images: "+"{:.3f}".format(t_filtering)+" seconds.",'blue','on_white',attrs=['bold']))
    if image_collection_filenames["SegPMRF"] is not None or image_collection["SegPMRF"] is not None:
        print(colored("Elapsed time to segment data using pmrf: "+"{:.3f}".format(t_pmrf)+" seconds.",'blue','on_white',attrs=['bold']))
    if image_collection_filenames["SegSRM"] is not None or image_collection["SegSRM"] is not None:
        print(colored("Elapsed time to segment data using srm: "+"{:.3f}".format(t_srm)+" seconds.",'blue','on_white',attrs=['bold']))
    if image_collection_filenames["SegKMEANS"] is not None or image_collection["SegKMEANS"] is not None:
        print(colored("Elapsed time to segment data using kmeans: "+"{:.3f}".format(t_kmeans)+" seconds.",'blue','on_white',attrs=['bold']))
    if image_collection_filenames["SegThresh"] is not None or image_collection["SegThresh"] is not None:
        print(colored("Elapsed time to segment data using threshold: "+"{:.3f}".format(t_thresh)+" seconds.",'blue','on_white',attrs=['bold']))
    if has_gt:
        print(colored("Elapsed time to calculate general metrics: "+"{:.3f}".format(t_general)+" seconds.",'blue','on_white',attrs=['bold']))
    if fiber_data:
        print(colored("Elapsed time to calculate fiber metrics: "+"{:.3f}".format(t_fibers)+" seconds.",'blue','on_white',attrs=['bold']))
    print(colored("Total elapsed time: "+"{:.3f}".format(t_total)+" seconds.",'blue','on_white',attrs=['bold']))
    print(colored("**********************************************************************************************************************************",'blue','on_white',attrs=['bold']))


if __name__ == "__main__":
    main()
