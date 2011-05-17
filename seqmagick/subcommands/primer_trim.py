"""
Find a primer sequence in a gapped alignment, trim to amplicon
"""
import argparse
import itertools

from Bio import SeqIO


def ungap_index(sequence, gap_chars='-'):
    """
    Iterates over a sequence, returning a the index for each position in the
    non-gapped sequence::

      raw_sequence[index] == raw_sequence.replace('-', '')[ungapped_index]

    Assuming index refers to a position:

    >>> ungap_index('AC-TG-')
    [0, 1, None, 2, 3, None]
    """
    counter = itertools.count(0)
    return [next(counter) if c not in gap_chars else None
            for c in iter(sequence)]


def build_parser(parser):
    parser.add_argument('source_file', help="Source alignment file",
            type=argparse.FileType('r'))
    parser.add_argument('output_file', help="Destination trimmed file",
            type=argparse.FileType('w'))
    parser.add_argument('forward_primer',
            help="The forward primer used")
    parser.add_argument('reverse_primer',
            help="The reverse primer used")
    parser.add_argument('--no-reverse-complement', default=False,
            action='store_true', help="""Do not take the reverse complement of
            the reverse primer""", dest="reverse_complement")
    parser.add_argument('--alignment-format', default='fasta',
            help='Alignment format (default: %(default)s)')


def action(arguments):
    """

    """
    # Load the alignment
    # AlignIO wo
    with arguments.source_file:
        sequences = SeqIO.parse(arguments.source_file,
                arguments.alignment_format)


    raise NotImplementedError()
