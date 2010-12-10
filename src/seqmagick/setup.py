#! /usr/bin/env python

import glob
from distutils.core import setup

setup(name = 'SeqMagick',
      version = '0.1',
      description = 'Tools for converting and working with sequence files.',
      author = 'Brian Hodges',
      author_email = 'bhodges@fhcrc.org',
      package_dir = {'seqmagick': '.'},
      packages = ['seqmagick'],
      scripts = glob.glob('scripts/*.py'),
      )
