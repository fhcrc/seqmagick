"""
Extract the sequence IDs from a file
"""
import argparse
import sys

from Bio import SeqIO

from seqmagick import fileformat

def build_parser(parser):
    parser.add_argument('sequence_file', help="Sequence file",
            type=argparse.FileType('r'))
    parser.add_argument('-o', '--output-file', help="Destination trimmed file",
            type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--source-format', default=None)

def action(arguments):
    # Determine file format for input and output
    source_format = (arguments.source_format or
            fileformat.from_filename(arguments.sequence_file.name))

    with arguments.sequence_file:
        ids = (sequence.id for sequence in SeqIO.parse(arguments.sequence_file,
                                                       source_format))
        with arguments.output_file:
            for i in ids:
                print >> arguments.output_file, i
