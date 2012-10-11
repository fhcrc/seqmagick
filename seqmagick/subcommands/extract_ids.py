"""
Extract the sequence IDs from a file
"""
import sys

from Bio import SeqIO

from seqmagick import fileformat

from . import common

def build_parser(parser):
    parser.add_argument('sequence_file', help="Sequence file",
            type=common.FileType('r'))
    parser.add_argument('-o', '--output-file', help="Destination trimmed file",
            type=common.FileType('w'), default=sys.stdout)
    parser.add_argument('--input-format', help="""Input format for sequence
            file""")
    parser.add_argument('-d', '--include-description', action='store_true',
            default=False, help="""Include the sequence description in output
            [default: %(default)s]""")

def action(arguments):
    common.exit_on_sigpipe()

    # Determine file format for input and output
    source_format = (arguments.input_format or
            fileformat.from_handle(arguments.sequence_file))

    with arguments.sequence_file:
        sequences = SeqIO.parse(arguments.sequence_file, source_format)
        if arguments.include_description:
            ids = (sequence.description for sequence in sequences)
        else:
            ids = (sequence.id for sequence in sequences)
        with arguments.output_file:
            for i in ids:
                print >> arguments.output_file, i
