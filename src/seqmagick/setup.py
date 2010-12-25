#! /usr/bin/env python

import glob
from distutils.core import setup

setup(name = 'seqmagick',
      version = '0.1',
      description = 'Tools for converting and working with sequence files.',
      author = 'Erick Matsen <ematsen@gmail.com>, Brian Hodges <bhodges@fhcrc.org>',
      package_dir = {'seqmagick': '.'},
      packages = ['seqmagick'],
      scripts = ['scripts/seqmagick'],
      )
      #scripts = glob.glob('scripts/*.py')'),
