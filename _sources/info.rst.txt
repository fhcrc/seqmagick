``info``
========

``seqmagick info`` describes one or more sequence files

Example
*******
::

    seqmagick info examples/*.fasta

    name                      alignment  min_len  max_len  avg_len  num_seqs
    examples/aligned.fasta    TRUE       9797     9797     9797.00  15
    examples/dewrapped.fasta  TRUE       240      240      240.00   148
    examples/range.fasta      TRUE       119      119      119.00   2
    examples/test.fasta       FALSE      972      9719     1573.67  15
    examples/wrapped.fasta    FALSE      120      237      178.50   2

Output can be in comma-separated, tab-separated, or aligned formats. See
``seqmagick info -h`` for details.

Usage:

.. literalinclude:: info.help
