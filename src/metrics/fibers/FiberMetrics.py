"""
.. module:: FiberMetrics
   :synopsis: This module is responsible for calculating the metric specific for fiber data. Metrics are calculated per slices and then per stack.

.. moduleauthor:: J. Mike Macneil and T. Perciano
"""

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "J. Mike Macneil and T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"

import matplotlib.pyplot as plt

from scipy import misc, ndimage as ndi
from scipy import interpolate
from skimage import io
from skimage import measure
from skimage.morphology import binary_erosion, convex_hull_image
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.util import img_as_uint, img_as_int, img_as_float
from skimage.measure import label, regionprops
from skimage.morphology import dilation, disk, erosion
from scipy.spatial import Voronoi, voronoi_plot_2d
from ...util.util import listdir_fullpath, create_dir
from scipy.stats import genextreme as gev
from colorama import init
from termcolor import colored

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, Image
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch

import os
from glob import glob

import numpy as np
import csv

def my_first_page(canvas,doc):
    """
    Boiler plate from reportlab documentation
    """
    PAGE_HEIGHT = defaultPageSize[1]
    PAGE_WIDTH = defaultPageSize[0]
    Title = "MSMcam Dashboard: Fiber Metrics Report"
    page_info = ""
    canvas.saveState()
    canvas.setFont('Times-Bold',16)
    canvas.drawCentredString(PAGE_WIDTH/2.0 , PAGE_HEIGHT-108, Title)
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75*inch, "First Page / %s" % page_info)
    canvas.restoreState()

def my_later_pages(canvas,doc):
    page_info = ""
    canvas.saveState()
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75 * inch, "Page %d %s" % (doc.page, page_info))
    canvas.restoreState()


def write_pdf_report(input_settings,run_number):
    styles = getSampleStyleSheet()
    
    input_type = input_settings["InputType"]
    input_dir = input_settings["InputDir"]
    if input_type==1:
        input_dir_split = input_dir.split("/")
        input_dir = "/".join(input_dir_split[:-1])+"/"
        input_file = input_dir
    
    output_dir = input_dir+"/msmcam_run_"+str(run_number)+'/metrics/dashboard/'
    figures_output_dir = input_dir+"/msmcam_run_"+str(run_number)+'/metrics/dashboard/figures/'

    PDFName = "Dashboard.pdf"

    doc = SimpleDocTemplate(output_dir+PDFName)
    Story = [Spacer(1,2*inch)]
    style = styles["Normal"]

    image_filenames = glob(figures_output_dir+"*.png")
    num_figs = len(image_filenames)
    print(num_figs-num_figs%2-1)

    for fig_number in range(num_figs-num_figs%2-1):
        plot_1 = Image(image_filenames[fig_number],width = 4*inch,height = 3.5*inch)
        plot_2 = Image(image_filenames[fig_number+1], width = 4*inch,height = 3.5*inch)
        Story.append(Table([[plot_1,plot_2]] ) )
        Story.append(Spacer(1,0.2*inch))

    doc.build(Story, onFirstPage=my_first_page, onLaterPages=my_later_pages)

def make_plots(input_settings, run_number):

    """
    Read in the fieldnames.
    """
    
    input_type = input_settings["InputType"]
    input_dir = input_settings["InputDir"]
    if input_type==1:
        input_dir_split = input_dir.split("/")
        input_dir = "/".join(input_dir_split[:-1])+"/"
        input_file = input_dir
    
    figures_output_dir = input_dir+"/msmcam_run_"+str(run_number)+'/metrics/dashboard/figures/'
    create_dir(figures_output_dir)
    alg_paths = glob(input_dir+"/msmcam_run_"+str(run_number)+'/metrics/fibers/*/')
    num_algs = len(alg_paths)
    names_algs = list()
    slice_metrics_data_matrix_list = list()

    for i in range(num_algs):
        if 'pmrf' in alg_paths[i]:
            names_algs.append('pmrf')
        if 'srm' in alg_paths[i]:
            names_algs.append('srm')
        if 'kmeans' in alg_paths[i]:
            names_algs.append('kmeans')
        if 'gt' in alg_paths[i]:
            names_algs.append('gt')
        if 'thresh' in alg_paths[i]:
            names_algs.append('thresh')
        fieldnames = []
        csv_filename = glob(alg_paths[i]+'SliceAnalysis.csv')[0]
        with open(csv_filename) as slicecsv:
            reader = csv.DictReader(slicecsv)
            fieldnames = reader.fieldnames
        """
        Reading in using loadtxt is just simpler than going row by row
        """
        slice_metrics_data_matrix_list.append(np.loadtxt(csv_filename, delimiter = ",",skiprows = 1))


    """
    Create a plot for each column
    """
    if len(slice_metrics_data_matrix_list[0])>2:
        for col_number, name in enumerate(fieldnames):
            plt.figure()
            plt.xlabel("Slice Number")
            plt.ylabel(name)
            plt.title("Slice Results for " + name)

            for i in range(len(names_algs)):
                data = np.asarray(slice_metrics_data_matrix_list[i][:,col_number])
                n_slices = len(data)
                x = np.asarray(range(n_slices))
                new_n_slices = np.linspace(x.min(),x.max(),n_slices*3)
                data_smooth = interpolate.spline(x,data,new_n_slices)
                plt.plot(new_n_slices, data_smooth)
            plt.legend(names_algs, loc='upper left')
            plt.savefig(figures_output_dir+name+".png")
            plt.close()

    return

class FiberMetrics:
    """ A class containing specific metrics for fiber data

    :param seg_image:       Segmented image
    :type seg_image:        ndimage
    :param input_settings:  Dictionary with input settings
    :type input_settings:   dic
    :param algorithm:       Algorithm name
    :type algorithm:        str

    """

    def __init__(self, seg_image, input_settings,algorithm,run_number):
        """
         Initialize FiberMetrics directories
        """
        self.seg_image = seg_image
        """ Segmented image """
        
        self.input_type = input_settings["InputType"]
        input_dir = input_settings["InputDir"]
        if self.input_type==1:
            input_dir_split = input_dir.split("/")
            input_dir = "/".join(input_dir_split[:-1])+"/"
            input_file = input_dir
        self.input_dir = input_dir
        
        self.output_dir = self.input_dir+"/msmcam_run_"+str(run_number)+'/metrics/fibers/'+algorithm
        """ Output directory """
        self.watershed_dir = os.path.join(self.output_dir,'watershed')
        """ Directory containing watershed results """
        self.voronoi_output_dir = os.path.join(self.output_dir,'voronoi')
        """ Directory containing voronoi results """
        self.voronoi_fib_output_dir = os.path.join(self.output_dir,'voronoi_fibers')
        """ Directory containing voronoi+fibers results overlap """
        self.slice_output_csv_filename = os.path.join(self.output_dir,'SliceAnalysis.csv')
        """ Name of the file with mtrics per slice """
        self.stack_output_csv_filename = os.path.join(self.output_dir,'StackAnalysis.csv')
        """ Name of the file with metrics per stack """
        self.in_memory = input_settings["InMemory"]
        
        create_dir(self.output_dir)
        create_dir(self.watershed_dir)
        create_dir(self.voronoi_output_dir)
        create_dir(self.voronoi_fib_output_dir)
        self.alg = algorithm
        """ Algorithm name """

    def segmentation_to_watershed_and_voronoi(self, img, local_max_min_distance = 6, local_voronoi_distance = 6):
        """
        Converts a binary segmentation into a watershedded
        image to break up globs. Amounts to porting
        the fiji script into python.
        Two major steps.
        1) Get the watershedded image.
        2) Get the mask, dilate the image then invert it.
        """

        """
        Create the Watershed binary
        """
        img_watershed = img.copy()

        distance_watershed = ndi.distance_transform_edt(img_watershed)
        local_maxi_watershed = peak_local_max(distance_watershed, min_distance = local_max_min_distance, indices = False, labels = img_watershed, exclude_border=False)
        markers_watershed = label(local_maxi_watershed)
        labels_watershed = watershed(-distance_watershed, markers_watershed, mask = img_watershed, watershed_line = True)
        output_watershed_binary = np.zeros_like(img,dtype = "uint8")
        output_watershed_binary[labels_watershed >= 1] = 255

        """
        Draw the Voronoi regions.
        """
        #: Now we want to separate the two objects in image
        #: Generate the markers as local maxima of the distance to the background
        img_voronoi = (255-img).copy()

        distance_voronoi = ndi.distance_transform_edt(img_voronoi)
        local_maxi_voronoi = peak_local_max(-distance_voronoi, min_distance = local_voronoi_distance, indices=False, exclude_border=False)
        labels_voro = watershed(distance_voronoi, labels_watershed, watershed_line = True)

        mask = self.get_mask(img)
        #plt.imshow(mask)
        #plt.show()

        return output_watershed_binary, labels_voro*mask


    def get_mask(self, img):
        """ Function to return mask for an image

        :param img: Image
        :type img:  ndimage
        :returns:   Mask for segmented image
        :rtype:     ndimage

        """
        mask = convex_hull_image(img)
        return mask


    def run_metrics(self, verbose = True):
        """
        Output:
         slice_metrics: a list of dictionaries,
         each dictionary containing the results for a given slice.
         stack_metrics: a dictionary with the aggregate results of the entire stack.


        Notes:
        1) Three main steps are open image,
           write out binary and fiber results, then compute metrics.
           Each step should be reasonably well marked with comments.
        2) Aggregate metrics are computed on the fly,
           rather than looping through at the end.
        3) The choice of saving as a list of dictionaries is mainly
           because that seems to be the easiest way to write to csv.

        TODO: Right now, list of metrics we're interested in is hard-coded, would be nice to
        figure out how to separate this from the implementation and processing.
        """

        slice_metrics = []
        stack_metrics = {}
        number_of_slices = None

        if self.in_memory:
            number_of_slices, y, x  = self.seg_image.shape
        else:
            y, x = io.imread(self.seg_image[0]).shape
            number_of_slices = len(self.seg_image)

        for slice_number in range(number_of_slices):
            slice_metrics.append({})

            if self.in_memory:
                img = self.seg_image[slice_number,:,:]
            else:
            	img = io.imread(self.seg_image[slice_number])
            img[img==1]=255
            """
            Get both Voronoi and Watershed
            """
            watershed_out, voronoi_out = self.segmentation_to_watershed_and_voronoi(img,8,8)

            temp = watershed_out.copy()
            index = np.where(temp==1)
            temp[index] = 255
            name = str(self.watershed_dir)+"/Slice"+"{0:0>4}".format(slice_number)+".tif"
            misc.imsave(name,temp)

            height, width = voronoi_out.shape
            """
            Convert voronoi input into binary to do region counting etc.
            """
            output_voronoi_filename = os.path.join(self.voronoi_output_dir,"4p3_Voronoi_"+str(slice_number+1).zfill(4)+".tif")
            voronoi_output = np.zeros_like(img,dtype=np.uint8)
            voronoi_output[voronoi_out>=1] = 255
            io.imsave(output_voronoi_filename,voronoi_output)
            """
            Extract fibers to publish fiber image.
            """
            voronoi_fibers = np.zeros((height,width),dtype = "uint8")
            voronoi_fibers[img==255] = 255
            #voronoi_fibers[watershed==1] = 255
            voronoi_fibers[voronoi_out==0] = 127

            voronoi_fibers_outputfilename = os.path.join(self.voronoi_fib_output_dir,"4p3_VF_"+ str(slice_number+1).zfill(4) +".tif")
            io.imsave(voronoi_fibers_outputfilename, voronoi_fibers)
            """
            Extract connected regions on top of the fiber regions.
            """
            region_properties = measure.regionprops(voronoi_out, watershed_out)
            number_of_regions = len(region_properties)
            if number_of_regions<3:
                return None, None
            """
            Metric Set One: Calculates the mean and std dev of the regions areas
            by looping through the connected components. These are used later for th
            Fraction of the unoccupied area as well. We figure out the number
            of occupied regions here to avoid looping twice
            """
            regions_area_all = np.zeros(number_of_regions-1)
            occupied = np.zeros(number_of_regions-1)

            for region_number, region in enumerate(region_properties):
                if region_number == 0:
                    continue
                else:
                    regions_area_all[region_number-1] = region.area
                    occupied[region_number-1] = region.mean_intensity*region.area

            slice_metrics[-1]["Regions Area Mean"] = np.mean(regions_area_all)
            slice_metrics[-1]["Regions Area Std Dev"] = np.std(regions_area_all)

            self.add_to_stack_metric("Stack Regions Area Mean",
                                slice_metrics[-1]["Regions Area Mean"],
                                number_of_slices,
                                stack_metrics)
            self.add_to_stack_metric("Stack Regions Area Std Dev",
                                slice_metrics[-1]["Regions Area Std Dev"],
                                number_of_slices,
                                stack_metrics)

            """
            Metric Set Two: The mean and std dev unoccupied fraction for each region.
            Uses the previously calculated occupied array and region_areas
            """
            frac_unocc_all = (regions_area_all - occupied)/regions_area_all
            slice_metrics[-1]["Fraction Unoccupied Mean"] = np.mean(frac_unocc_all)
            slice_metrics[-1]["Fraction Unoccupied Std Dev"] = np.std(frac_unocc_all)

            self.add_to_stack_metric("Stack Fraction Unoccupied Mean",
                                slice_metrics[-1]["Fraction Unoccupied Mean"],
                                number_of_slices,
                                stack_metrics)
            self.add_to_stack_metric("Stack Fraction Unoccupied Std Dev",
                                slice_metrics[-1]["Fraction Unoccupied Std Dev"],
                                number_of_slices,
                                stack_metrics)

            """
            Metric Set Three: The mean and standard deviation
            of area weighted fraction unoccupied for each cell.
            """
            area_unocc_all = regions_area_all*frac_unocc_all
            slice_metrics[-1]["Area Unoccupied Mean"] = np.mean(area_unocc_all)
            slice_metrics[-1]["Area Unoccupied Std Dev"] = np.std(area_unocc_all)

            self.add_to_stack_metric("Stack Area Unoccupied Mean",
                                slice_metrics[-1]["Area Unoccupied Mean"],
                                number_of_slices,
                                stack_metrics)
            self.add_to_stack_metric("Stack Area Unoccupied Std Dev",
                                slice_metrics[-1]["Area Unoccupied Std Dev"],
                                number_of_slices,
                                stack_metrics)
            """
            Metric Set Four: The mean and standard deviation of the area
            weighted normalized unoccupied fraction for each
            region.
            """
            total_area = np.sum(regions_area_all)
            area_unocc_norm_all = (regions_area_all*frac_unocc_all/total_area)
            slice_metrics[-1]["Area Unoccupied Normalized Mean"] = np.mean(area_unocc_norm_all)
            slice_metrics[-1]["Area Unoccupied Normalized Std Dev"] = np.std(area_unocc_norm_all)

            self.add_to_stack_metric("Stack Area Unoccupied Normalized Mean",
                                slice_metrics[-1]["Area Unoccupied Normalized Mean"],
                                number_of_slices,
                                stack_metrics)
            self.add_to_stack_metric("Stack Area Unoccupied Normalized Std Dev",
                                slice_metrics[-1]["Area Unoccupied Normalized Std Dev"],
                                number_of_slices,
                                stack_metrics)
            """
            Metric Set Five: Alpha Parameter and Alpha unoccupied parameter.
            """
            slice_metrics[-1]["Alpha"] = (np.mean(regions_area_all*regions_area_all)
                                          /(np.mean(regions_area_all)**2))
            slice_metrics[-1]["Alpha Unoccupied"] = (np.mean(area_unocc_all*area_unocc_all)
                                                     /(np.mean(area_unocc_all)**2))

            self.add_to_stack_metric("Stack Alpha",
                                slice_metrics[-1]["Alpha"],
                                number_of_slices,
                                stack_metrics)
            self.add_to_stack_metric("Stack Alpha Unoccupied",
                                slice_metrics[-1]["Alpha Unoccupied"],
                                number_of_slices,
                                stack_metrics)
            """
            Metric Set Six: Calculates the porosity
            """
            temp_frac_unocc = frac_unocc_all.copy()
            temp_frac_unocc[temp_frac_unocc >= 0.98] = float("nan")

            ones = np.ones_like(temp_frac_unocc)

            slice_metrics[-1]["Mean Porosity"] = (np.nansum(regions_area_all*temp_frac_unocc)
                                                  / total_area)
            slice_metrics[-1]["Mean Porosity for F"] = (np.nansum(regions_area_all*temp_frac_unocc
                                                          /(ones - temp_frac_unocc)) / total_area)
            slice_metrics[-1]["Mean Porosity for K0"] = (np.nansum(regions_area_all*temp_frac_unocc**2
                                                           /((ones-temp_frac_unocc)**2)) / total_area)

            self.add_to_stack_metric("Stack Mean Porosity",
                                slice_metrics[-1]["Mean Porosity"],
                                number_of_slices,
                                stack_metrics)
            self.add_to_stack_metric("Stack Mean Porosity for F",
                                slice_metrics[-1]["Mean Porosity for F"],
                                number_of_slices,
                                stack_metrics)
            self.add_to_stack_metric("Stack Mean Porosity for K0",
                                slice_metrics[-1]["Mean Porosity for K0"],
                                number_of_slices,
                                stack_metrics)

            if verbose:
                print(colored('Completed metrics for slice number '+str(slice_number+1),
                              'yellow','on_grey', attrs=['bold']))

        """
        End of Main Loop
        """

        """
        Polish off the rest of the stack metrics.
        """
        self.complete_stack_metrics(slice_metrics,stack_metrics)

        return slice_metrics, stack_metrics

    def complete_stack_metrics(self,slice_metrics, stack_metrics):
        """
        A couple of the metrics computed in the paper aren't
        very simple to do on the fly, so we have this final helper function
        to polish it off.

        :param slice_metrics:   Dictionary with slice metrics
        :type slice_metrics:    dic
        :param stack_metrics:   Dictionary with stack metrics
        :type stack_metrics:    dic

        """

        stack_metrics["Stack Regions Area Fractional Std Dev"] = (stack_metrics["Stack Regions Area Std Dev"]
                                                                  /stack_metrics["Stack Regions Area Mean"])
        stack_metrics["Stack Fraction Unoccupied Fractional Std Dev"] = (stack_metrics["Stack Fraction Unoccupied Std Dev"]
                                                                         /stack_metrics["Stack Fraction Unoccupied Mean"])
        stack_metrics["Stack Area Unoccupied Fractional Std Dev"] = (stack_metrics["Stack Area Unoccupied Std Dev"]
                                                                     /stack_metrics["Stack Area Unoccupied Mean"])
        stack_metrics["Stack Area Unoccupied Normalized Fractional Std Dev"] = (stack_metrics["Stack Area Unoccupied Std Dev"]
                                                                               /stack_metrics["Stack Area Unoccupied Mean"])
        number_of_slices = len(slice_metrics)

        mean_porosity = np.zeros(number_of_slices)
        mean_porosity_for_F = np.zeros(number_of_slices)
        mean_porosity_for_K0 = np.zeros(number_of_slices)

        for slice_data in slice_metrics:
            mean_porosity = slice_data["Mean Porosity"]
            mean_porosity_for_F = slice_data["Mean Porosity for F"]
            mean_porosity_for_K0 = slice_data["Mean Porosity for K0"]


        stack_metrics["Stack Mean Porosity Std Dev"] = np.std(mean_porosity)
        stack_metrics["Stack Mean Porosity for F Std Dev"] = np.std(mean_porosity)
        stack_metrics["Stack Mean Porosity for K0 Std Dev"] = np.std(mean_porosity)

        return

    def add_to_stack_metric(self,metric_key, slice_metric, number_of_slices, stack_metrics):
        """
        Simple helper to add to the stack metrics, since the if and else statements
        clutter up the main loop.

        """

        if not metric_key in stack_metrics.keys():
            stack_metrics[metric_key] = slice_metric/number_of_slices
        elif metric_key in stack_metrics.keys():
            stack_metrics[metric_key] += slice_metric/number_of_slices

        return

    def publish_results(self,slice_metrics,stack_metrics):
        """
        This is a little naughty, but currently hardcoding the metric fields.
        This should be kept in mind in case I want to update this.
        Split the writing into two steps. One for slice metrics and
        another for stack metrics.
        """
        with open(self.slice_output_csv_filename, "w") as slicecsv:
            fieldnames = ["Regions Area Mean",
                          "Regions Area Std Dev",
                          "Fraction Unoccupied Mean",
                          "Fraction Unoccupied Std Dev",
                          "Area Unoccupied Mean",
                          "Area Unoccupied Std Dev",
                          "Area Unoccupied Normalized Mean",
                          "Area Unoccupied Normalized Std Dev",
                          "Alpha",
                          "Alpha Unoccupied",
                          "Mean Porosity",
                          "Mean Porosity for F",
                          "Mean Porosity for K0"]

            writer = csv.DictWriter(slicecsv,fieldnames = fieldnames)

            writer.writeheader()
            for line in slice_metrics:
                writer.writerow(line)


        with open(self.stack_output_csv_filename, "w") as stackcsv:
            fieldnames = ["Stack Regions Area Mean",
                          "Stack Regions Area Std Dev",
                          "Stack Regions Area Fractional Std Dev",
                          "Stack Fraction Unoccupied Mean",
                          "Stack Fraction Unoccupied Std Dev",
                          "Stack Fraction Unoccupied Fractional Std Dev",
                          "Stack Area Unoccupied Mean",
                          "Stack Area Unoccupied Std Dev",
                          "Stack Area Unoccupied Fractional Std Dev",
                          "Stack Area Unoccupied Normalized Mean",
                          "Stack Area Unoccupied Normalized Std Dev",
                          "Stack Area Unoccupied Normalized Fractional Std Dev",
                          "Stack Alpha",
                          "Stack Alpha Unoccupied",
                          "Stack Mean Porosity",
                          "Stack Mean Porosity Std Dev",
                          "Stack Mean Porosity for F",
                          "Stack Mean Porosity for F Std Dev",
                          "Stack Mean Porosity for K0",
                          "Stack Mean Porosity for K0 Std Dev"]

            writer = csv.DictWriter(stackcsv,fieldnames = fieldnames)

            writer.writeheader()
            writer.writerow(stack_metrics)

        return

    def process(self):
        """ Calculates metrics for all the segmentation algorithms """
        if self.alg=='pmrf':
            print(colored('PMRF','yellow','on_grey', attrs=['bold']))
        if self.alg=='srm':
            print(colored('SRM','yellow','on_grey', attrs=['bold']))
        if self.alg=='kmeans':
            print(colored('KMEANS','yellow','on_grey', attrs=['bold']))
        if self.alg=='gt':
            print(colored('GT','yellow','on_grey', attrs=['bold']))
        if self.alg=='thresh':
            print(colored('Threshold','yellow','on_grey', attrs=['bold']))
        slice_metrics, stack_metrics = self.run_metrics()
        if slice_metrics is None:
            return False
        self.publish_results(slice_metrics, stack_metrics)
        return True
