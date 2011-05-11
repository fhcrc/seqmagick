#! /usr/bin/env python

from distutils.core import setup

setup(name = 'seqmagick',
      version = '0.1',
      description = 'Tools for converting and working with sequence files.',
      author = 'Erick Matsen <ematsen@gmail.com>, Brian Hodges <bhodges@fhcrc.org>',
      packages = ['seqmagick', 'seqmagick.scripts', 'seqmagick.test'],
      scripts = ['scripts/seqmagick'],
      )
