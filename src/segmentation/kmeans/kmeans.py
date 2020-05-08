"""
    .. module:: Kmeans
    :synopsis: This module is responsible for image segmentation using kmeans algorithm
    
    .. moduleauthor:: H. Krishnan
    
    """

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "H. Krishnan"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"


from ...util.util import create_dir, upscale, listdir_fullpath
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin
from sklearn.datasets import load_sample_image
from sklearn.utils import shuffle
from time import time
from scipy import misc
from skimage import io, measure
from colorama import init
from termcolor import colored
from ..SegmentationAlgorithm import SegmentationAlgorithm
import os

class KMEANS(SegmentationAlgorithm):

    def __init__(self, image, input_settings, preproc_settings,seg_settings,run_number):
        super().__init__(image, input_settings, preproc_settings,seg_settings,run_number)

    def recreate_image(self, codebook, labels, w, h):
        """Recreate the (compressed) image from the code book & labels
            
        :param codebook:        Clusters centers
        :type codebook:         array
        :param labels:          Labels for each pixel
        :type labels:           array
        :param w:               Image width
        :type w:                int
        :param h:               Image height
        :type h:                int
        :returns:               Recreated classified image
        :rtype:                 ndimage

        """
        d = codebook.shape[1]
        image = np.zeros((w, h, d))
        image = np.zeros((w, h))
        label_idx = 0
        for i in range(w):
            for j in range(h):
                image[i][j] = codebook[labels[label_idx]]
                label_idx += 1
        #print image.shape
        return image

    def segment(self):
	    """ Main segmentation process using kmeans algorithm 
	        
	    """
	
	    input_type = self.input_settings["InputType"]
	    input_dir = self.input_settings["InputDir"]
	    if input_type==1:
	        input_dir_split = input_dir.split("/")
	        input_dir = "/".join(input_dir_split[:-1])
	        input_file = input_dir
	
	    num_labels = self.seg_settings["NumClustersKMEANS"]
	    multiphase = self.seg_settings["Multiphase"]
	    masked = self.input_settings["Masked"]
	    invert = self.seg_settings["Invert"]
	    in_memory = self.input_settings["InMemory"]
	    first_slice = self.input_settings["FirstSlice"]
	    last_slice = self.input_settings["LastSlice"]
	    num_slices = last_slice-first_slice
	    if in_memory:
	        self.segmented = np.zeros(self.image.shape,dtype=np.uint8)
	
	    if multiphase is False:
	        if masked:
	            n_colors = 3
	        else:
	            n_colors = 2
	    else:
	        n_colors = num_labels

	    # create outout directories, if writing tifs
	    if not in_memory:
	        output_tif_dir = os.path.join(input_dir,"msmcam_run_"+str(self.run_number),"res/kmeans/tif")
	        create_dir(output_tif_dir)

	    #if in_memory:
	    #    self.image = upscale(self.image, self.preproc_settings["DownsampleScale"])
	    j = first_slice
	    for i in range(num_slices):
	        if in_memory:
	            img = self.image[i,:,:]
	        else:
	            img = io.imread(str(self.image[i]))
	        
	        image_array = img
	        img = np.array(img, dtype=np.float64) / 255

	        w, h = original_shape = tuple(img.shape)
	        d = 1
	        image_array = np.reshape(image_array, (w * h, d))

	        #print("Fitting model on a small sub-sample of the data")
	        t0 = time()
	        image_array_sample = shuffle(image_array, random_state=0)[:1000]
	        kmeans = KMeans(n_clusters=n_colors, random_state=0).fit(image_array_sample)
	        #print("done in %0.3fs." % (time() - t0))

	        # Get labels for all points
	        #print("Predicting color indices on the full image (k-means)")
	        t0 = time()
	        labels = kmeans.predict(image_array)
	        #print("done in %0.3fs." % (time() - t0))

	        rec_img = self.recreate_image(kmeans.cluster_centers_, labels, w, h)
	        
	        if multiphase is False:
	            res = labels.astype(np.uint8)
	            res.shape = img.shape
	            if masked is False:
	                index1 = np.where(res==0)
	                index2 = np.where(res==1)
	            else:
	                label_background = res[0,0]
	                print(label_background)
	                count = 0
	                for label in range(n_colors):
	                    if label != label_background:
	                        if count==0:
	                            index1 = np.where(res==label)
	                        if count==1:
	                            index2 = np.where(res==label)
	                        count = count+1
	                res[res==label_background]=0
	            if len(index1[0])>len(index2[0]):
	                res[index1]=255
	                res[index2]=0
	                if invert:
	                    res[index1]=0
	                    res[index2]=255
	            else:
	                res[index2]=255
	                res[index1]=0
	                if invert:
	                    res[index2]=0
	                    res[index1]=255

	            if self.input_settings['FiberData']:
	                concomp = measure.label(res)
	                props = measure.regionprops(concomp)
	                for k in range(len(props)):
	                    if props[k].eccentricity>0.97 or props[k].area<100:
	                        label = props[k].label
	                        res[concomp==label]=0
	            
	            if in_memory:
	                self.segmented[i,:,:] = res
	            else:
	                name = str(output_tif_dir)+"/Slice"+"{0:0>4}".format(j)+".tif"
	                misc.imsave(name,res)
	                self.segmented_filenames.append(name)
	        else:
	            if in_memory:
	                self.segmented[i,:,:] = rec_img
	            else:
	                name = str(output_tif_dir)+"/Slice"+"{0:0>4}".format(j)+".tif"
	                misc.imsave(name,rec_img)
	                self.segmented_filenames.append(name)
	        print(colored('Finished image '+str(j),'yellow','on_grey', attrs=['bold']))
	        j = j+1

