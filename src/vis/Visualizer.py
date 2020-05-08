"""
    .. module:: Visualizer
    :synopsis: This module is responsible for providing visualization tools for MSMcam

    .. moduleauthor:: D. Ushizima and T. Perciano

    """

__copyright__   =   "CAMERA Materials Segmentation & Metrics (MSM) Copyright (c) 2017, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved."
__developers__  =   "D. Ushizima, H. Krishnan, T. Perciano, M. MacNeil, D. Parkinson"
__author__      =   "D. Ushizima and T. Perciano"
__maintainer__  =   "CAMERA Image Processing"
__email__       =   "camera.image.processing@gmail.com"
__license__     =   "Modified BSD"
__version__     =   "0.1"


import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')
from skimage import color, io, img_as_ubyte, external
from skimage.transform import resize
from skimage.filters import threshold_otsu, sobel
from skimage.morphology import closing, square, remove_small_objects, thin
from ..util.util import create_dir, listdir_fullpath
import imageio
from colorama import init
from termcolor import colored
import numpy as np
import os

class Visualizer:
    """ A Visualizer

    :param input_dir:       Input directory
    :type input_dir:        str
    :param vis_settings:    Dictionary with visualization parameters
    :type vis_settings:     dic

    """

    def __init__(self, input_dir, input_type, vis_settings, in_memory, run_number):
        """ Initialize Visualizer

        """
        self.input_type = input_type
        self.input_dir = input_dir
        """ Input directory """
        if self.input_type==1:
            input_dir_split = input_dir.split("/")
            self.input_dir = "/".join(input_dir_split[:-1])+"/"
            self.input_file = input_dir
            """ Input file names if image stack """
        self.vis_settings = vis_settings
        """ Dictionary with visualization settings """
       	self.in_memory = in_memory
       	self.run_number = run_number

    def viewGeneralMetrics(self, df, output_dir):
        """ Plots graphics with summary for general metrics

        :param df:          Dataframe with general metrics
        :type df:           pandas dataframe
        :param output_dir:  Output directory
        :type output_dir:   str

        """
        #fig, axes = plt.subplots(nrows=2,ncols=2,figsize=(16,6))
        fig = plt.figure(figsize=(16,9))
        plt.suptitle('Summary of accuracy results',fontsize=15)

        ax1 = plt.subplot2grid((2, 2), (0, 0))
        df[['Precision','Recall','Accuracy','Specificity','F','SegPorosity', 'GtPorosity']].plot(ax=ax1,kind='bar')
        plt.sca(ax1)
        plt.xticks(rotation='horizontal')
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax1.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05),fancybox=True, shadow=True,ncol=4)

        ax2 = plt.subplot2grid((2, 2), (0, 1))
        df[['SegVolume','GtVolume']].plot(ax=ax2,kind='bar')
        plt.sca(ax2)
        plt.xticks(rotation='horizontal')
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05),fancybox=True, shadow=True,ncol=4)

        ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=2)
        ax3.axis("off")
        table = plt.table(cellText=df.values, rowLabels = df.index, colLabels = df.columns,loc='center')
        table.set_fontsize(11)
        table.scale(1.1, 1.1)
        if self.vis_settings["SavePlots"]:
            fig.savefig(output_dir+'summaryplot.png')
        if self.vis_settings["ViewPlots"]:
            plt.show()

    def viewSegmentations(self, image_collection, image_collection_filenames):
        """ View overlay with segmentation results

        :param image_collection:    File names for all the segmentation results
        :type image_collection:     list

        """
        minSize = 50
        out_path = self.input_dir+"msmcam_run_"+str(self.run_number)+"/vis"
        if self.in_memory is False:
            create_dir(out_path)
        #TODO: (2) what if the processed slices are smaller than slices available
        ic_list = list()

        overlay_filenames = list()
        overlay = None

        if self.in_memory:
            ic_orig = image_collection["Orig"]
            ic_list.append(ic_orig)
            ic_gt = image_collection["GT"]
            if ic_gt is not None:
                ic_list.append(ic_gt)
            ic_res_pmrf = image_collection["SegPMRF"]
            if ic_res_pmrf is not None:
                ic_list.append(ic_res_pmrf)
            ic_res_kmeans = image_collection["SegKMEANS"]
            if ic_res_kmeans is not None:
                ic_list.append(ic_res_kmeans)
            ic_res_srm = image_collection["SegSRM"]
            if ic_res_srm is not None:
                ic_list.append(ic_res_srm)
            ic_res_thresh = image_collection["SegThresh"]
            if ic_res_thresh is not None:
                ic_list.append(ic_res_thresh)
            z, h, w = ic_orig.shape
            overlay = np.zeros((h,w,3,z), dtype=np.uint8)
        else:
            ic_orig = image_collection_filenames["Orig"]
            ic_list.append(ic_orig)
            ic_gt = image_collection_filenames["GT"]
            if ic_gt is not None:
                ic_list.append(ic_gt)
            ic_res_pmrf = image_collection_filenames["SegPMRF"]
            if ic_res_pmrf is not None:
                ic_list.append(ic_res_pmrf)
            ic_res_kmeans = image_collection_filenames["SegKMEANS"]
            if ic_res_kmeans is not None:
                ic_list.append(ic_res_kmeans)
            ic_res_srm = image_collection_filenames["SegSRM"]
            if ic_res_srm is not None:
                ic_list.append(ic_res_srm)
            ic_res_thresh = image_collection_filenames["SegThresh"]
            if ic_res_thresh is not None:
                ic_list.append(ic_res_thresh)
            im = io.imread(str(ic_orig[0]))
            h, w = im.shape
            z = len(image_collection_filenames["Filtered"])
        fh = h/1000 #anim size"
        if h<1000:
            fh = 1
        x=round(h/fh)
        y=round(w/fh)

        #new_list = [ len(w) for w in ic_list]
        #minLen = min(new_list)

        #for each image, for each method, paint the edges accordingly
        for n in range(z):
            if self.in_memory:
                im_orig = ic_orig[n,:,:]
            else:
                im_orig = io.imread(str(ic_orig[n]))
            im_orig = img_as_ubyte(im_orig)
            im_orig= resize(im_orig,(x,y),mode='reflect')
            im_orig = color.grey2rgb(im_orig, alpha=None)

            if ic_gt is not None:
                if self.in_memory:
                    im = img_as_ubyte(ic_gt[n,:,:])
                else:
                    im = img_as_ubyte(io.imread(str(ic_gt[n])))
                im = resize(im,(x,y),mode='reflect')
                try:
                    thresh = threshold_otsu(im)  #needed due to interpolation problems
                    im = closing(im > thresh, square(3))
                    im = remove_small_objects(im, minSize)
                except:
                    pass
                edg = sobel(im)
                edg = edg >0
                edg = thin(edg)
                im_orig[edg, 0] = 1    # Red
                #if darker borders
                im_orig[edg, 1] = 0
                im_orig[edg, 2] = 0

            if ic_res_kmeans is not None:
                if self.in_memory:
                    im = img_as_ubyte(ic_res_kmeans[n,:,:])
                else:
                    im = img_as_ubyte(io.imread(str(ic_res_kmeans[n])))
                im = resize(im,(x,y),mode='reflect')
                try:
                    thresh = threshold_otsu(im)  #needed due to interpolation problems
                    im = closing(im > thresh, square(3))
                    im = remove_small_objects(im, minSize)
                except:
                    pass
                edg = sobel(im)
                edg = edg >0
                edg = thin(edg)
                im_orig[edg, 1] = 1     # Green
                #if darker borders
                im_orig[edg, 0] = 0
                im_orig[edg, 2] = 0

            if ic_res_srm is not None:
                if self.in_memory:
                    im = img_as_ubyte(ic_res_srm[n,:,:])
                else:
                    im = img_as_ubyte(io.imread(str(ic_res_srm[n])))
                im = resize(im,(x,y),mode='reflect')
                try:
                    thresh = threshold_otsu(im)
                    im = closing(im > thresh, square(3))
                    im = remove_small_objects(im, minSize)
                except:
                    pass
                edg = sobel(im)
                edg = edg >0
                edg = thin(edg)
                im_orig[edg, 2] = 1     # Blue
                im_orig[edg, 0] = 0
                im_orig[edg, 1] = 0

            if ic_res_pmrf is not None:
                if self.in_memory:
                    im = img_as_ubyte(ic_res_pmrf[n,:,:])
                else:
                    im = img_as_ubyte(io.imread(str(ic_res_pmrf[n])))
                im = resize(im,(x,y),mode='reflect')
                try:
                    thresh = threshold_otsu(im)
                    im = closing(im > thresh, square(3))
                    im = remove_small_objects(im, minSize)
                except:
                    pass
                edg = sobel(im)
                edg = edg >0
                edg = thin(edg)
                im_orig[edg, 0] = 1     # Dark orange
                im_orig[edg, 1] = 0.5
                im_orig[edg, 2] = 0
                
            if ic_res_thresh is not None:
                if self.in_memory:
                    im = img_as_ubyte(ic_res_thresh[n,:,:])
                else:
                    im = img_as_ubyte(io.imread(str(ic_res_thresh[n])))
                im = resize(im,(x,y),mode='reflect')
                edg = sobel(im)
                edg = edg >0
                edg = thin(edg)
                im_orig[edg, 0] = 0.75     # Lilac
                im_orig[edg, 1] = 0.6
                im_orig[edg, 2] = 0.75

            im_orig = img_as_ubyte(im_orig)
            if self.in_memory:
                overlay[:,:,:,n] = im_orig
            else:
                #saving slices
                name = str(out_path)+"/Slice"+"{0:0>4}".format(n)+".tif"
                external.tifffile.imsave(name, im_orig,photometric='rgb')
                overlay_filenames.append(name)
                
            print(colored('Finished image '+str(n),'yellow','on_grey', attrs=['bold']))
        return overlay, overlay_filenames

    def create_gif(self):
        print(colored('Creating .gif animation...','green','on_grey', attrs=['bold']))
        filenames = sorted(listdir_fullpath(str(self.input_dir+"msmcam_run_"+str(self.run_number)+"/vis/")))
        if '.DS_Store' in filenames[0]:
            filenames = filenames[1:]
        movie_path = self.input_dir+"msmcam_run_"+str(self.run_number)+ "/vis/movie"
        create_dir(movie_path)
        with imageio.get_writer(os.path.join(movie_path,'movie.gif'), mode='I',duration=0.5) as writer:
            for filename in filenames:
                image = imageio.imread(filename)
                writer.append_data(image)
