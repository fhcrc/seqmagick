#! /usr/bin/env python

from distutils.core import setup
import sys

from seqmagick import __version__ as version

if sys.version_info < (2, 7):
    print 'ERROR: seqmagick requires at least Python 2.7 to run.'
    sys.exit(1)

try:
    import Bio
except:
    print "ERROR: BioPython isn't installed."
    sys.exit(1)

setup(name='seqmagick',
      version=version,
      description='Tools for converting and modifying sequence files '
                  'from the command-line',
      url='http://github.com/fhcrc/seqmagick',
      author='Matsen Group',
      author_email='http://matsen.fhcrc.org/',
      #author_email="matsen@fhcrc.org",
      packages=['seqmagick', 'seqmagick.scripts', 'seqmagick.test',
                'seqmagick.subcommands'],
      scripts=['scripts/seqmagick'],
      requires=['Python (>= 2.7)', 'biopython (>=1.54)'],
      classifiers=[
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
          ],
      license="GPL V3",
      )
