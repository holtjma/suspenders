from setuptools import setup
import MergeImprove

setup(name='suspenders',
      version='0.1.2',
      description='Allows the merging of alignments that have been annotated using pylapels into a single alignment that picks the highest quality alignment.',
      url='http://code.google.com/p/suspenders',
      author='James Holt',
      author_email='holtjma@cs.unc.edu',
      license='MIT',
      packages=['MergeImprove'],
      install_requires=['pysam'],
      scripts=['bin/pysuspenders'],
      zip_safe=False)
