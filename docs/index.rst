.. seqmagick documentation master file, created by
   sphinx-quickstart on Thu May 19 16:18:13 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. Contents:


.. "Fork me on github"


.. raw:: html

    <a href="http://github.com/fhcrc/seqmagick"><img style="position: absolute; top: 0; left: 0; border: 0;" src="_static/fork.png" alt="Fork me on GitHub"></a>


=========
seqmagick
=========

.. contents::
   :depth: 4
   :class: new

.. toctree::
   :maxdepth: 1

   changelog


Motivation
==========

We often have to convert between sequence formats and do little tasks on them,
and it's not worth writing scripts for that.  Seqmagick is a kickass little
utility built in the spirit of imagemagick_ to expose the file format
conversion in Biopython in a convenient way.  Instead of having a big mess of
scripts, there is one that takes arguments::

    seqmagick convert a.fasta b.phy    # convert from fasta to phylip
    seqmagick mogrify --ungap a.fasta  # remove all gaps from a.fasta, in place
    seqmagick info *.fasta             # describe all FASTA files in the current directory

And more.

Installation
============

First, you'll need to install `BioPython`_. NumPy (which parts of BioPython
depend on) is not required for ``seqmagick`` to function. Once done, install
the latest release with::

    pip install seqmagick

Or install the bleeding edge version::

    pip install git+git://github.com/fhcrc/seqmagick.git@master#egg-info=seqmagick

Use
===

Seqmagick can be used to query information about sequence files, convert
between types, and modify sequence files.  All functions are accessed through
subcommands::

    seqmagick <subcommand> [options] arguments

List of Subcommands
===================

.. toctree::
   :maxdepth: 2

   convert_mogrify
   backtrans_align
   extract_ids
   info
   quality_filter
   primer_trim

Supported File Extensions
=========================

By default, ``seqmagick`` infers the file type from extension. Currently mapped
extensions are:

.. include:: extensions.rst

.. note::

    NEXUS-format output requires the ``--alphabet`` flag.

Default Format
--------------

When reading from stdin or writing to stdout, ``seqmagick`` defaults to fasta
format.  This behavior may be overridden with the ``--input-format`` and
``--output-format`` flags.

If an extension is not listed, you can either rename the file to a supported
extension, or specify it manually via ``--input-format`` or ``--output-format``.

Compressed file support
-----------------------

most commands support gzip (files ending in ``.gz``) and bzip (files ending in
``.bz2`` or ``.bz``) compressed inputs and outputs. File types for these files
are inferred using the extension of the file after stripping the file extension
indicating that the file is compressed, so ``input.fasta.gz`` would be inferred
to be in FASTA format.

Acknowledgements
================

seqmagick is written and maintained by the `Matsen Group`_ at the Fred
Hutchinson Cancer Research Center.


Contributing
============

We welcome contributions! Simply fork the repository `on GitHub`_ and send a pull request.

.. _`on GitHub`: http://github.com/fhcrc/seqmagick/
.. _`Matsen Group`: http://matsen.fhcrc.org/
.. _imagemagick: http://www.imagemagick.org/script/command-line-tools.php
.. _`BioPython`: http://www.biopython.org/
