``protparam``
===========================

``protparam`` calculates the molecular weight (Mw) and theoretical
isoelectric point (pI) for a set of protein sequences.

Examples
--------

Basic calculation
^^^^^^^^^^^^^^^^^

Take two sequences like `haemoglobins.fasta`_::

    $ seqmagick protparam haemoglobins.fasta

Gives the tab-delimited output::

    sp|P69905|HBA_HUMAN	15256.89	8.72
    sp|P68871|HBB_HUMAN	15997.81	6.74

where the first column is the sequence id, the second the is calculated
molecular weight in daltons and the third is the predicted isoelectric
point (pI).

Sorting
^^^^^^^

Output can also be sorted based on molecular weight, isoelectric point,
sequence length or alphabetically by name using the ``--sort`` option.

To sort based on molecular weight, descending::

    $ seqmagick protparam --sort mass-desc haemoglobins.fasta
    sp|P68871|HBB_HUMAN	15997.81	6.74
    sp|P69905|HBA_HUMAN	15256.89	8.72

.. note::
  ``protparam`` expects ungapped sequences with uppercase letters.
  If your sequences contain gaps or lowercase letters, you should
  fix this using ``convert`` or ``mogrify`` with the ``--ungapped``
  or ``--upper`` options, respectively.

  Only the 20 standard amino acids are supported by default.
  If your sequence contains characters for unknown amino acids (eg
  'X', 'B'), you can use the ``--allow-unknown-residues`` to treat
  them as untitratable with respect to pI calculation, with an average
  mass of 112.5 Da.


Support for other physiochemical properties computed by the
`BioPython ProtParam module`_ or it's web-based inspiration
`ExPASy ProtParam`_ are not yet implemented.

Command-line Arguments
**********************

.. program-output:: ../seqmagick.py protparam -h

.. _`haemoglobins.fasta`: http://www.uniprot.org/uniprot/?query=id:P68871+OR+id:P69905&format=fasta

.. _`BioPython ProtParam module`: http://biopython.org/wiki/ProtParam

.. _`ExPASy ProtParam`: http://web.expasy.org/protparam/