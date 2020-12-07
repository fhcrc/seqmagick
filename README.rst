=========
seqmagick
=========

.. image:: https://travis-ci.org/fhcrc/seqmagick.svg?branch=master
    :target: https://travis-ci.org/fhcrc/seqmagick

We often have to convert sequence files between formats and do little
manipulations on them, and it's not worth writing scripts for that.
``seqmagick`` is a kickass little utility to expose the file format conversion
in BioPython in a convenient way.  Instead of having a big mess of scripts,
there is one that takes arguments::

    seqmagick convert a.fasta b.phy    # convert from fasta to phylip
    seqmagick mogrify --ungap a.fasta  # remove all gaps from a.fasta, in place
    seqmagick info *.{fasta,sto}       # describe all FASTA and Stockholm
                                       # files in the current directory

Requirements
============

* Python >= 3.5
* biopython >= 1.78

Installation
============

Use pip::

   pip install seqmagick

Note that as of version 0.8.0, this package requires Python 3.5+. If
you want to use the most recent version compatible with Python 2.7::

  pip install seqmagick==0.6.2

Features
========

* Modifying sequences: Remove gaps, reverse complement, reverse, change case,

  - Remove gaps
  - Reverse & reverse complement
  - Trim to a range of residues
  - Change case
  - Sort by length or ID
  - `more`_

* Displaying `information <http://seqmagick.readthedocs.org/en/latest/info.html>`_ about
  sequence files
* Subsetting sequence files by:

  - Position
  - ID
  - Deduplication
  - `more`_

* Filtering sequences by `quality score
  <http://seqmagick.readthedocs.org/en/latest/quality_filter.html>`_
* Trimming alignments to a `region of interest
  <http://seqmagick.readthedocs.org/en/latest/primer_trim.html>`_ defined by the
  forward and reverse primers

Want to learn more? Head to the `Documentation`_.

``seqmagick`` is free software under the GPL v3.


.. _`Documentation`: http://seqmagick.readthedocs.org/en/latest/

.. _`more`: http://seqmagick.readthedocs.org/en/latest/convert_mogrify.html
