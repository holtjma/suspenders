from setuptools import setup

setup(name='suspenders',
      version='0.1',
      description='Allows the merging of alignments that have been annotated using pylapels into a single alignment that picks the highest quality alignment.',
      url='<INSERT URL>',
      author='James Holt',
      author_email='holtjma@cs.unc.edu',
      license='MIT',
      packages=['MergeImprove'],
      install_requires=['pysam'],
      scripts=['bin/pysuspenders'],
      zip_safe=False)