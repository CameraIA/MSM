

from scipy import ndimage as ndi
import numpy as np

from skimage import io
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from skimage.util import img_as_uint, img_as_int
from skimage.measure import label, regionprops
from skimage.morphology import dilation, disk, erosion


"""
Module of functions to be able to bring segmentation results
into 
"""


def segmentation_to_watershed_and_voronoi(img,
                                          local_max_min_distance = 4,
                                          local_voronoi_distance = 4):
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
    local_maxi_watershed = peak_local_max(distance_watershed,
                                min_distance = local_max_min_distance,
                                indices = False,
                                labels = img_watershed)
    markers_watershed = label(local_maxi_watershed)
    labels_watershed = watershed(-distance_watershed,
                       markers_watershed,
                       mask = img_watershed,
                       watershed_line = True)
    output_watershed_binary = np.zeros_like(img,dtype = "uint8")
    output_watershed_binary[labels_watershed >= 1] = 1

    region_properties_watershed = regionprops(labels_watershed)

    """
    Draw the Voronoi regions.
    """
    # Now we want to separate the two objects in image
    # Generate the markers as local maxima of the distance to the background
    img_voronoi = (255-img).copy()

    distance_voronoi = ndi.distance_transform_edt(img_voronoi)
    local_maxi_voronoi = peak_local_max(-distance_voronoi, min_distance =
                                        local_voronoi_distance, indices=False)
    labels_voro = watershed(distance_voronoi, labels_watershed, watershed_line = True)

    mask = get_mask(img)

    return output_watershed_binary, labels_voro*mask
    

def get_mask(img, dilation_disk_size = 5):
    """
    Simple method to extract a mask from the segmentation.
    main idea is to dilate the daylights out of the segmentation and 
    then take the inverse.
    """

    dilated = dilation(img,disk(dilation_disk_size))
    labelled_inverted = label(1-dilated)

    output = np.zeros_like(dilated,dtype = "uint8")
    output[labelled_inverted != 1] = 1
   
    
    return output



