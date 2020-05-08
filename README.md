# MSM
Materials Segmentation & Metrics: tools to process image segmentation using traditional codes and evaluation metrics to check uncertainty

MSM package as described in the Journal of Synchrotron Radiation article:
"Insight into 3D micro-CT data: exploring segmentation algorithms through performance metrics"
by
Perciano,  Ushizima,  Krishnan,  Parkinson,  Larson,  Pelt, Bethel, Zok and Sethian, although the codes also count with contributions from J. Mike Macneil.

[Full paper] http://onlinelibrary.wiley.com/doi/10.1107/S1600577517010955/abstract

### What is this repository for? ###

* MSMcam
* 0.1

### How do I start? ###

* Install: **sudo python setup.py install** or **python setup.py install --user**
* Change the following variables within MSMSettings.ini:
    * InputDir					=	/yourfolderwith/msm/examples/data/gambier_region/orig/
    * GroundTruthDir				=	/yourfolderpath/msm/examples/data/gambier_region/groundtruth/
* Run: MSMcam -c MSMSettings.ini

### How do I process images from a different folder?

* Update MSMSettings.ini as before with the path to your original images as well as groundtruth

### Python versions we tested on:
* 2.7
* 3.4
* 3.5

### General notes:

* We suggest converting pixel depth to 8-bit, but MSM will take care of this conversion if needed;

### Known limitations:
* Multiphase samples with more than 2 phases might present innacuracies;
* PMRF function is not available yet;

### FAQ:
1) I could not install using sudo python setup.py install. What else should I try?
A = It might be the case that you need to install the packages manually using conda or pip. This might solve known problems during installation on linux systems through easy_install

