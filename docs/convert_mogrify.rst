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

    seqmagick mogrify --tail 5 --reverse-complement --ungap --deduplicate-sequences examples/test.fasta examples/test.fasta

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

The full set of options to ``mogrify`` and ``convert`` are:

Sequence File Modification
^^^^^^^^^^^^^^^^^^^^^^^^^^
::

      --line-wrap N         Adjust line wrap for sequence strings. When N is 0,
                            all line breaks are removed. Only fasta files are
                            supported for the output format.
      --sort {length-asc,length-desc,name-asc,name-desc}
                            Perform sorting by length or name, ascending or
                            descending. ASCII sorting is performed for names

Sequence Modification
^^^^^^^^^^^^^^^^^^^^^
::

      --apply-function /path/to/module.py:function_name
                            Specify a custom function to apply to the input
                            sequences, specified as
                            /path/to/file.py:function_name. Function should accept
                            an iterable of Bio.SeqRecord objects, and yield
                            SeqRecords. Specify more than one to chain.
      --cut start:end[,start2:end2]
                            1-indexed start and end positions for cutting sequences, : separated. Includes last item. Start or end can be left unspecified
                            to indicate start/end of sequence.
      --relative-to ID      Apply --cut relative to the indexes of non-gap residues in sequence identified by ID
      --dash-gap            Change . and : into - for all sequences
      --mask start:end[,start2:end2...]
                        Replace residues in 1-indexed slice with gap-
                        characters. If --relative-to is also specified,
                        coordinates are relative to the sequence ID provided.
      --lower               Translate the sequences to lower case
      --reverse             Reverse the order of sites in sequences
      --reverse-complement  Convert sequences into reverse complements
      --squeeze             Remove any gaps that are present in the same position
                            across all sequences in an alignment (equivalent to
                            --squeeze-threshold=1.0)
      --squeeze-threshold PROP
                            Trim columns from an alignment which have gaps in
                            least the specified proportion of sequences.
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

Record Selection
^^^^^^^^^^^^^^^^
::

      --deduplicate-sequences
                            Remove any duplicate sequences by sequence content,
                            keep the first instance seen
      --deduplicated-sequences-file FILE
                            Write all of the deduplicated sequences to a file
      --deduplicate-taxa    Remove any duplicate sequences by ID, keep the first
                            instance seen
      --exclude-from-file FILE
                            Filter sequences, removing those sequence IDs in the
                            specified file
      --include-from-file FILE
                            Filter sequences, keeping only those sequence IDs in
                            the specified file
      --head N              Trim down to top N sequences
      --max-length N        Discard any sequences beyond the specified maximum
                            length. This operation occurs *before* all length-
                            changing options such as cut and squeeze.
      --min-length N        Discard any sequences less than the specified minimum
                            length. This operation occurs *before* all length-
                            changing options such as cut and squeeze.
      --min-ungapped-length N
                            Discard any sequences less than the specified minimum
                            length, excluding gaps. This operation occurs *before*
                            all length-changing options such as cut and squeeze.
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

Sequence ID Modification
^^^^^^^^^^^^^^^^^^^^^^^^
::

      --first-name          Take only the first whitespace-delimited word as the
                            name of the sequence
      --name-suffix SUFFIX  Append a suffix to all IDs.
      --name-prefix PREFIX  Insert a prefix for all IDs.
      --pattern-replace search_pattern replace_pattern
                            Replace regex pattern "search_pattern" with
                            "replace_pattern" in sequence ID
      --strip-range         Strip ranges from sequences IDs, matching </x-y>

Format Options
^^^^^^^^^^^^^^

By default, file format is inferred from extension::


      --input-format Format
                            Input file format (default: determine from extension)
      --output-format Format
                            Output file format (default: determine from extension)


.. _`BioPython SeqIO wiki page`: http://www.biopython.org/wiki/SeqIO#File_Formats
