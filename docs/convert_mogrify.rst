``convert`` and ``mogrify``
===========================

Convert and mogrify achieve similar goals. ``convert`` performs some operation
on a file (from changing format to something more complicated) and writes to a
new file. ``mogrify`` modifies a file in place, and would not normally be used
to convert formats.

The two have similar signatures::

    seqmagick convert [options] infile outfile

vs::

    seqmagick mogrify [options] infile

Options are shared between convert and mogrify.

Examples
--------

Basic Conversion
^^^^^^^^^^^^^^^^

``convert`` can be used to convert between any file types BioPython supports
(which is many). For a full list of supported types, see the `BioPython SeqIO
wiki page`_.

By default, file type is inferred from file extension, so::

    seqmagick convert a.fasta a.sto

converts an existing file ``a.fasta`` from FASTA to Stockholm format. **Neat!**
But there's more.

Sequence Modification
^^^^^^^^^^^^^^^^^^^^^

A wealth of options await you when you're ready to do something slightly more
complicated with your sequences.

Let's say I just want a few of my sequences::

    $ seqmagick convert --head 5 examples/test.fasta examples/test.head.fasta
    $ seqmagick info examples/test*.fasta
    name                      alignment  min_len  max_len  avg_len  num_seqs
    examples/test.fasta       FALSE      972      9719     1573.67  15
    examples/test.head.fasta  FALSE      978      990      984.00   5

Or I want to remove any gaps, reverse complement, select the last 5 sequences,
and remove any duplicates from an alignment in place::

    seqmagick mogrify --tail 5 --reverse-complement --ungap --deduplicate-sequences examples/test.fasta

You can even define your own functions in python and use them via
``--apply-function``.

.. note::
  To maximize flexibility, most transformations passed as options to
  ``mogrify`` and ``convert`` are processed *in order*, so::

       seqmagick convert --min-length 50 --cut 1:5 a.fasta b.fasta

  will work fine, but::

       seqmagick convert --cut 1:5 --min-length 50 a.fasta b.fasta

  will never return records, since the cutting transformation happens before
  the minimum length predicate is applied.

Command-line Arguments
**********************

.. literalinclude:: convert.help

.. _`BioPython SeqIO wiki page`: http://www.biopython.org/wiki/SeqIO#File_Formats
