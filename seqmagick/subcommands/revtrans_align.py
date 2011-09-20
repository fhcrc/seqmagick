"""
Given a protein alignment and FASTA file, align the FASTA file
"""

# TODO: Add tests
# TODO: Specify input and output gap characters
# TODO: Infer output format from extension, default to fasta

import argparse
import itertools
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqmagick import fileformat

from . import common

def build_parser(parser):
    parser.add_argument('protein_align', help="""Protein Alignment""")
    parser.add_argument('nucl_align', help="""FASTA Alignment""")
    parser.add_argument('--out-file', type=argparse.FileType('w'),
            default=sys.stdout, metavar='destination_file', help="""Output
            destination. Default: STDOUT""")

    return parser

def batch(iterable, size):
    i = iter(iterable)
    while True:
        yield list(itertools.islice(i, size))

def convert(prot_seq, nucl_seq):
    # Ungap nucleotides
    codons = batch(str(nucl_seq.seq.ungap('-')), 3)
    codons = (''.join(i) for i in codons)
    # TODO: Check for valid conversion in translation table
    nucl_align = ['---' if p == '-' else next(codons)
                  for p in str(prot_seq.seq)]
    try:
        next(codons)
        raise ValueError("Additional codons present")
    except:
        # OK
        pass
    return SeqRecord(Seq(''.join(nucl_align)), id=nucl_seq.id,
            description=nucl_seq.description)

def convert_all(prot_alignment, nucl_sequences):
    """
    Convert protein sequences to nucleotide alignment
    """
    zipped = itertools.izip_longest(prot_alignment, nucl_sequences)
    for p, n in zipped:
        if p is None:
            raise ValueError("Exhausted protein sequences")
        elif n is None:
            raise ValueError("Exhausted nucleotide sequences")
        yield convert(p, n)

def action(arguments):
    """
    Run
    """
    # Ignore SIGPIPE, for head support
    common.exit_on_sigpipe()

    prot_sequences = SeqIO.parse(arguments.protein_align,
            fileformat.from_filename(arguments.protein_align))
    nucl_sequences = SeqIO.parse(arguments.nucl_align,
            fileformat.from_filename(arguments.nucl_align))

    SeqIO.write(convert_all(prot_sequences, nucl_sequences),
            arguments.out_file, 'fasta')
