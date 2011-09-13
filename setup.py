#! /usr/bin/env python

# Use setuptools, falling back on provide
try:
    from setuptools import setup, find_packages
except ImportError:
    import distribute_setup
    distribute_setup.use_setuptools()
    from setuptools import setup, find_packages

import sys
from distutils import log

from seqmagick import __version__ as version

if sys.version_info < (2, 7):
    print 'ERROR: seqmagick requires at least Python 2.7 to run.'
    sys.exit(1)

requires = ['biopython>=1.58']

setup(name='seqmagick',
      version=version,
      description='Tools for converting and modifying sequence files '
                  'from the command-line',
      url='http://github.com/fhcrc/seqmagick',
      download_url='http://pypi.python.org/pypi/seqmagick',
      author='Matsen Group',
      author_email='http://matsen.fhcrc.org/',
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'seqmagick = seqmagick.scripts.cli:main'
      ]},
      install_requires=requires,
      classifiers=[
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      license="GPL V3",
      )
