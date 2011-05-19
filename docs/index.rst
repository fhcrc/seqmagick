.. seqmagick documentation master file, created by
   sphinx-quickstart on Thu May 19 16:18:13 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========
seqmagick
=========

Contents:

.. toctree::
   :maxdepth: 2

Motivation
==========

We often have to convert between sequence formats and do little tasks on them,
and it's not worth writing scripts for that.  Seqmagick is a kickass little
utility built in the spirit of imagemagick_ to expose the file format
conversion in Biopython in a convenient way.  Instead of having a big mess of
scripts, there is one that takes arguments.

Installation
============

First, you'll need to install `BioPython`_. NumPy (which parts of BioPython
depend on) is not required for ``seqmagick`` to function. Once done, install
with::

    pip install seqmagick

Get the bleeding edge version `here
<https://github.com/fhcrc/seqmagick/zipball/master>`_, or clone our
repository::

    git clone git://github.com/fhcrc/seqmagick.git

Use
===

Seqmagick can be used to query information about sequence files, convert
between types, and modify sequence files.  All functions are accessed through
``seqmagick <subcommand>``

Subcommands
===========

``info``
========

``seqmagick info`` describes one or more sequence files::

    $ seqmagick info examples/*.fasta

    name                      alignment  min_len  max_len  avg_len  num_seqs
    examples/aligned.fasta    TRUE       9797     9797     9797.00  15
    examples/dewrapped.fasta  TRUE       240      240      240.00   148
    examples/range.fasta      TRUE       119      119      119.00   2
    examples/test.fasta       FALSE      972      9719     1573.67  15
    examples/wrapped.fasta    FALSE      120      237      178.50   2

``convert`` and ``mogrify``
===========================

Convert and mogrify achieve similar goals. ``convert`` performs some operation
on a file (from changing format to something more complicated) and writes to a
new file. ``mogrify`` modifies a file in place, and would not normally be used
to convert formats.

Basic Conversion
----------------

seqmagick can be used to convert between any two file types BioPython
supported. For a full list of supported types, see the `BioPython SeqIO wiki page`_.

By default, file type is inferred from file extension, so::

    seqmagick convert a.fasta a.sto

converts an existing file ``a.fasta`` from FASTA to Stockholm format. Neat! But there's more.

Sequence Modification
---------------------

A wealth of options await you when you're ready to do something slightly more complicated with your sequences.

Let's say I want a few of my sequences::

    $ seqmagick convert --head 5 examples/test.fasta examples/test.head.fasta
    $ seqmagick info examples/test*.fasta
    name                      alignment  min_len  max_len  avg_len  num_seqs
    examples/test.fasta       FALSE      972      9719     1573.67  15
    examples/test.head.fasta  FALSE      978      990      984.00   5

Or I want to remove any gaps, reverse complement, select the last 5 sequences, and remove any duplicates from an alignment in place::

    seqmagick mogrify --tail 5 --reverse-complement --ungap --deduplicate-sequences examples/test.fasta examples/test.fasta

You can even define your own functions in python and use them via ``--apply-function``.

The full set of options to ``mogrify`` and ``convert`` are:

Sequence File Modification::

      --line-wrap N         Adjust line wrap for sequence strings. When N is 0,
                            all line breaks are removed. Only fasta files are
                            supported for the output format.
      --sort {length-asc,length-desc,name-asc,name-desc}
                            Perform sorting by length or name, ascending or
                            descending. ASCII sorting is performed for names

Sequence Modificaton::

      --apply-function /path/to/module.py:function_name
                            Specify a custom function to apply to the input
                            sequences, specified as
                            /path/to/file.py:function_name. Function should accept
                            an iterable of Bio.SeqRecord objects, and yield
                            SeqRecords. Specify more than one to chain.
      --cut start:end       1-indexed start and end positions for cutting
                            sequences, : separated. Includes last item.
      --dash-gap            Change . and : into - for all sequences
      --lower               Translate the sequences to lower case
      --reverse             Reverse the order of sites in sequences
      --reverse-complement  Convert sequences into reverse complements
      --squeeze             Remove any gaps that are present in the same position
                            across all sequences in an alignment
      --transcribe {dna2rna,rna2dna}
                            Transcription and back transcription for generic DNA
                            and RNA. Source sequences must be the correct alphabet
                            or this action will likely produce incorrect results.
      --translate {dna2protein,rna2protein,dna2proteinstop,rna2proteinstop}
                            Translate from generic DNA/RNA to proteins. Options
                            with "stop" suffix will NOT translate through stop
                            codons .Source sequences must be the correct alphabet
                            or this action will likely produce incorrect results.
      --ungap               Remove gaps in the sequence alignment
      --upper               Translate the sequences to upper case

Record Selection::

      --deduplicate-sequences
                            Remove any duplicate sequences by sequence content,
                            keep the first instance seen
      --deduplicated-sequences-file FILE
                            Write all of the deduplicated sequences to a file
      --deduplicate-taxa    Remove any duplicate sequences by ID, keep the first
                            instance seen
      --head N              Trim down to top N sequences
      --max-length N        Discard any sequences beyond the specified maximum
                            length. This operation occurs *before* all length-
                            changing options such as cut and squeeze.
      --min-length N        Discard any sequences less than the specified minimum
                            length. This operation occurs *before* all length-
                            changing options such as cut and squeeze.
      --pattern-include regex
                            Filter the sequences by regular expression in name
      --pattern-exclude regex
                            Filter out sequences by regular expression in name
      --prune-empty         Prune sequences containing only gaps ('-')
      --seq-pattern-include regex
                            Filter the sequences by regular expression in sequence
      --seq-pattern-exclude regex
                            Filter out sequences by regular expression in sequence
      --tail N              Trim down to bottom N sequences

Sequence ID Modification::

      --first-name          Take only the first whitespace-delimited word as the
                            name of the sequence
      --name-suffix SUFFIX  Append a suffix to all IDs.
      --name-prefix PREFIX  Insert a prefix for all IDs.
      --pattern-replace search_pattern replace_pattern
                            Replace regex pattern "search_pattern" with
                            "replace_pattern" in sequence ID
      --strip-range         Strip ranges from sequences IDs, matching </x-y>

Format Options::

      --input-format Format
                            Input file format (default: determine from extension)
      --output-format Format
                            Output file format (default: determine from extension)


``primer-trim``
===============

``primer-trim`` trims an alignment to a region defined by a set of forward and
reverse primers.

Possibly to implement:
----------------------
check
  check integrity of files
  should take --DNA --AA arguments to check if the file is one of these types
  (or perhaps an alphabet)
  make sure that you are using IUPAC characters for AA or Nuc


.. _imagemagick: http://www.imagemagick.org/script/command-line-tools.php

.. _`BioPython SeqIO wiki page`: http://www.biopython.org/wiki/SeqIO#File_Formats
.. _`BioPython`: http://www.biopython.org/


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

