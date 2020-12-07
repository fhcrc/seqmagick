#!/usr/bin/env python

import os
import sys
import subprocess
from setuptools import setup, find_packages

subprocess.call(
    ('mkdir -p seqmagick/data && '
     'git describe --tags --dirty > seqmagick/data/ver.tmp '
     '&& mv seqmagick/data/ver.tmp seqmagick/data/ver '
     '|| rm -f seqmagick/data/ver.tmp'),
    shell=True, stderr=open(os.devnull, "w"))

# must import __version__ after call to 'git describe' above
from seqmagick import __version__

setup(name='seqmagick',
      version=__version__,
      description='Tools for converting and modifying sequence files '
      'from the command-line',
      url='http://github.com/fhcrc/seqmagick',
      download_url='http://pypi.python.org/pypi/seqmagick',
      author='Matsen Group',
      # author_email='http://matsen.fhcrc.org/',
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'seqmagick = seqmagick.scripts.cli:main'
          ]},
      package_data={
          'seqmagick': ['data/*'],
          'seqmagick.test.integration': ['data/*']
      },
      setup_requires=['nose>=1.0'],
      python_requires='>=3.5',
      test_suite='nose.collector',
      install_requires=['biopython>=1.78'],
      classifiers=[
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Development Status :: 4 - Beta',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      license="GPL V3")
