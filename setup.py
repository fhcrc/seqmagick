#! /usr/bin/env python

from distutils.core import setup

from seqmagick import __version__ as version

setup(name='seqmagick',
      version=version,
      description='Tools for converting and working with sequence files.',
      author='Erick Matsen <ematsen@gmail.com>, Brian Hodges <bhodges@fhcrc.org>',
      packages=['seqmagick', 'seqmagick.scripts', 'seqmagick.test',
                'seqmagick.subcommands'],
      scripts=['scripts/seqmagick'],
      )
