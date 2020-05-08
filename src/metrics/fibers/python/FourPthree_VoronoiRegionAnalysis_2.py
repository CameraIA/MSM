
""" 
Python port of Natalies Voronoi Region Analysis Matlab Script.
"""

import numpy

from skimage import io
from skimage import measure
from skimage.util import img_as_float

from scipy.stats import genextreme as gev

import os
from glob import glob
import configparser
import csv

from SegmentationToVoronoi import segmentation_to_watershed_and_voronoi



def run_metrics(settings, verbose = True):
    """
    Parameters:
     settings: a dictionary indicating where the files
     containing the voronoi images and watershed results are.
     
     verbose: outputting to the command line progress. True by default.
      In the future should probably be better as a logger.

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

    Todo: Right now, list of metrics we're interested in is hard-coded, would be nice to
    figure out how to separate this from the implementation and processing.
    """

    slice_metrics = []
    stack_metrics = {}

    image_filenames =  glob(settings["input path"]+"*.tif")


    number_of_slices =len(image_filenames)

    for slice_number in range(number_of_slices):
        slice_metrics.append({})

        img = io.imread(image_filenames[slice_number])

        """
        Get both Voronoi and Watershed
        """
        watershed, voronoi = segmentation_to_watershed_and_voronoi(img)

        height, width = voronoi.shape
        """
        Convert voronoi input into binary to do region counting etc.
        """
        output_voronoi_filename = (settings["voronoi output path"]
                                   +"4p3_Voronoi_"+"%04d"%(slice_number+1)+".tif")
        voronoi_output = numpy.zeros_like(img)
        voronoi_output[voronoi>=1] = 255
        io.imsave(output_voronoi_filename,voronoi_output)

        """
        Extract fibers to publish fiber image.
        """
        voronoi_fibers = numpy.zeros((height,width),dtype = "uint8")
        for i in range(height):
            for j in range(width):
                if watershed[i,j] == 1:
                    voronoi_fibers[i,j] = 255
                if voronoi[i,j] == 0:
                    voronoi_fibers[i,j] = 127

        voronoi_fibers_outputfilename = (settings["voronoi fibers output path"] + "4p3_VF_"
                                         + "%04d"%(slice_number+1)+".tif")
        io.imsave(voronoi_fibers_outputfilename, voronoi_fibers)

        """
        Extract connected regions on top of the fiber regions.
        """
        region_properties = measure.regionprops(voronoi, watershed)
        number_of_regions = len(region_properties)

        """
        Metric Set One: Calculates the mean and std dev of the regions areas
        by looping through the connected components. These are used later for th 
        Fraction of the unoccupied area as well. We figure out the number
        of occupied regions here to avoid looping twice
        """
        regions_area_all = numpy.zeros(number_of_regions-1)
        occupied = numpy.zeros(number_of_regions-1)

        for region_number, region in enumerate(region_properties):
            if region_number == 0:
                continue
            else:
                regions_area_all[region_number-1] = region.area
                occupied[region_number-1] = region.mean_intensity*region.area

        slice_metrics[-1]["Regions Area Mean"] = numpy.mean(regions_area_all)
        slice_metrics[-1]["Regions Area Std Dev"] = numpy.std(regions_area_all)

        add_to_stack_metric("Stack Regions Area Mean",
                            slice_metrics[-1]["Regions Area Mean"],
                            number_of_slices,
                            stack_metrics)
        add_to_stack_metric("Stack Regions Area Std Dev",
                            slice_metrics[-1]["Regions Area Std Dev"],
                            number_of_slices,
                            stack_metrics)

        """
        Metric Set Two: The mean and std dev unoccupied fraction for each region.
        Uses the previously calculated occupied array and region_areas
        """
        frac_unocc_all = (regions_area_all - occupied)/regions_area_all
        slice_metrics[-1]["Fraction Unoccupied Mean"] = numpy.mean(frac_unocc_all)
        slice_metrics[-1]["Fraction Unoccupied Std Dev"] = numpy.std(frac_unocc_all)

        add_to_stack_metric("Stack Fraction Unoccupied Mean",
                            slice_metrics[-1]["Fraction Unoccupied Mean"],
                            number_of_slices,
                            stack_metrics)
        add_to_stack_metric("Stack Fraction Unoccupied Std Dev",
                            slice_metrics[-1]["Fraction Unoccupied Std Dev"],
                            number_of_slices,
                            stack_metrics)

        """
        Metric Set Three: The mean and standard deviation 
        of area weighted fraction unoccupied for each cell.
        """
        area_unocc_all = regions_area_all*frac_unocc_all
        slice_metrics[-1]["Area Unoccupied Mean"] = numpy.mean(area_unocc_all) 
        slice_metrics[-1]["Area Unoccupied Std Dev"] = numpy.std(area_unocc_all)

        add_to_stack_metric("Stack Area Unoccupied Mean",
                            slice_metrics[-1]["Area Unoccupied Mean"],
                            number_of_slices,
                            stack_metrics)
        add_to_stack_metric("Stack Area Unoccupied Std Dev",
                            slice_metrics[-1]["Area Unoccupied Std Dev"],
                            number_of_slices,
                            stack_metrics)        
        """
        Metric Set Four: The mean and standard deviation of the area 
        weighted normalized unoccupied fraction for each
        region.
        """
        total_area = numpy.sum(regions_area_all)
        area_unocc_norm_all = (regions_area_all*frac_unocc_all/total_area)
        slice_metrics[-1]["Area Unoccupied Normalized Mean"] = numpy.mean(area_unocc_norm_all)
        slice_metrics[-1]["Area Unoccupied Normalized Std Dev"] = numpy.std(area_unocc_norm_all)

        add_to_stack_metric("Stack Area Unoccupied Normalized Mean",
                            slice_metrics[-1]["Area Unoccupied Normalized Mean"],
                            number_of_slices,
                            stack_metrics)
        add_to_stack_metric("Stack Area Unoccupied Normalized Std Dev",
                            slice_metrics[-1]["Area Unoccupied Normalized Std Dev"],
                            number_of_slices,
                            stack_metrics)    

        """
        Metric Set Five: Alpha Parameter and Alpha unoccupied parameter.
        """
        slice_metrics[-1]["Alpha"] = (numpy.mean(regions_area_all*regions_area_all)
                                      /(numpy.mean(regions_area_all)**2))
        slice_metrics[-1]["Alpha Unoccupied"] = (numpy.mean(area_unocc_all*area_unocc_all)
                                                 /(numpy.mean(area_unocc_all)**2))

        add_to_stack_metric("Stack Alpha",
                            slice_metrics[-1]["Alpha"],
                            number_of_slices,
                            stack_metrics)
        add_to_stack_metric("Stack Alpha Unoccupied",
                            slice_metrics[-1]["Alpha Unoccupied"],
                            number_of_slices,
                            stack_metrics) 

        """
        Metric Set Six: Calculates the porosity 
        """
        temp_frac_unocc = frac_unocc_all.copy()
        temp_frac_unocc[temp_frac_unocc >= 0.98] = float("nan")

        ones = numpy.ones_like(temp_frac_unocc)

        slice_metrics[-1]["Mean Porosity"] = (numpy.nansum(regions_area_all*temp_frac_unocc)
                                              / total_area)
        slice_metrics[-1]["Mean Porosity for F"] = (numpy.nansum(regions_area_all*temp_frac_unocc
                                                      /(ones - temp_frac_unocc)) / total_area)
        slice_metrics[-1]["Mean Porosity for K0"] = (numpy.nansum(regions_area_all*temp_frac_unocc**2
                                                       /((ones-temp_frac_unocc)**2)) / total_area)

        add_to_stack_metric("Stack Mean Porosity",
                            slice_metrics[-1]["Mean Porosity"],
                            number_of_slices,
                            stack_metrics)
        add_to_stack_metric("Stack Mean Porosity for F",
                            slice_metrics[-1]["Mean Porosity for F"],
                            number_of_slices,
                            stack_metrics)
        add_to_stack_metric("Stack Mean Porosity for K0",
                            slice_metrics[-1]["Mean Porosity for K0"],
                            number_of_slices,
                            stack_metrics)


        print("Completed metrics for slice number "+str(slice_number+1))

    """
    End of Main Loop
    """

    """
    Polish off the rest of the stack metrics.
    """
    complete_stack_metrics(slice_metrics, stack_metrics)

    return slice_metrics, stack_metrics


def complete_stack_metrics(slice_metrics, stack_metrics):
    """
    A couple of the metrics computed in the paper aren't
    very simple to do on the fly, so we have this final helper function
    to polish it off.
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

    mean_porosity = numpy.zeros(number_of_slices)
    mean_porosity_for_F = numpy.zeros(number_of_slices)
    mean_porosity_for_K0 = numpy.zeros(number_of_slices)
    
    for slice_data in slice_metrics:
        mean_porosity = slice_data["Mean Porosity"]
        mean_porosity_for_F = slice_data["Mean Porosity for F"]
        mean_porosity_for_K0 = slice_data["Mean Porosity for K0"]


    stack_metrics["Stack Mean Porosity Std Dev"] = numpy.std(mean_porosity)
    stack_metrics["Stack Mean Porosity for F Std Dev"] = numpy.std(mean_porosity)
    stack_metrics["Stack Mean Porosity for K0 Std Dev"] = numpy.std(mean_porosity)

    return

def add_to_stack_metric(metric_key, slice_metric, number_of_slices, stack_metrics):
    """
    Simple helper to add to the stack metrics, since the if and else statements
    clutter up the main loop.
    """

    if not metric_key in stack_metrics.keys():
        stack_metrics[metric_key] = slice_metric/number_of_slices
    elif metric_key in stack_metrics.keys():
        stack_metrics[metric_key] += slice_metric/number_of_slices

    return

def publish_results(slice_metrics,stack_metrics, settings):
    """
    This is a little naughty, but currently hardcoding the metric fields. 
    This should be kept in mind in case I want to update this.
    Split the writing into two steps. One for slice metrics and
    another for stack metrics.
    """

    slice_filename = settings["slice output csv filename"]
    stack_filename = settings["stack output csv filename"]
    
    
    with open(slice_filename, "w") as slicecsv:
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


    with open(stack_filename, "w") as stackcsv:
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

def load_settings():
    """
    Helper function to parse the INI file.
    Read in directories as a dictionary
    so we can access files and folders directly.
    Does a quick check to see if output paths exist,
    and if not, creates them.
    """
    config = configparser.ConfigParser()
    config.read("MetricSettings.ini")

    settings = {}

    for key in config["Paths"]:
        settings[key] = config["Paths"].get(key)
        if "output" in key and "path" in key:
            if not os.path.exists(settings[key]):
                os.makedirs(settings[key])


    """
    Quickly make the dirs for output if they don't exist.
    """    
    return settings
    
if __name__ == "__main__":
    """
    Main execution point of the program
    """
    settings = load_settings()
    slice_metrics, stack_metrics = run_metrics(settings)
    publish_results(slice_metrics, stack_metrics, settings)
