import subprocess
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from distutils.core import setup, Extension

class CustomDevelopCommand(develop):
    """
    Custom develop setup creating/removing a symlink to src.
    This is necessary because of this unsolved problem of setuptools:
    https://bitbucket.org/pypa/setuptools/issue/230
    """
    def run(self):
        if self.uninstall:
            self.multi_version = True
            self.uninstall_link()
            try:
                subprocess.check_output(["rm", "MSMcam"], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                pass
        else:
            self.install_for_development()
            self.warn_deprecated_options()
            try:
                subprocess.check_output(["ln", "-s", "src", "MSMcam"], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError:
                pass

setup(name='MSMcam',
      version='0.1',
      description='',
      url='git clone git@bitbucket.org:quantct/msm.git',
      author='CAMERA CIP',
      author_email='camera.image.processing@gmail.com',
      license='Modified BSD',
      cmdclass={'develop':CustomDevelopCommand},
      package_dir={'MSMcam':'src'},
      package_data={'MSMcam':['segmentation/srm/pysrm/*.so', 'segmentation/srm/pysrm/*.dylib', 'segmentation/pmrf/*.so', 'segmentation/pmrf/*.dylib']},
      packages=['MSMcam',
      			'MSMcam.segmentation',
                'MSMcam.segmentation.kmeans',
                'MSMcam.segmentation.pmrf',
                'MSMcam.segmentation.srm',
                'MSMcam.segmentation.srm.pysrm',
                'MSMcam.segmentation.threshold',
                'MSMcam.metrics',
                'MSMcam.metrics.fibers',
                'MSMcam.metrics.general',
                'MSMcam.inout',
                'MSMcam.util',
                'MSMcam.preprocessing',
                'MSMcam.postprocessing',
                'MSMcam.vis'],
      entry_points = {
          'console_scripts':
          ['MSMcam=MSMcam.MSMcam:main',]
          },
      install_requires=["numpy>=1.13.0", "scipy==0.19.1", "pandas", "xlsxwriter", "scikit-image==0.13.0", "scikit_learn>=0.19", "termcolor", "reportlab==3.4.0", "imageio", "matplotlib>=1.5.3"],
      zip_safe=False)
