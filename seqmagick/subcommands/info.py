"""
Info action
"""

import argparse
import collections
import csv
import os.path
import sys

from Bio import SeqIO

from seqmagick import fileformat


def build_parser(parser):
    parser.add_argument('source_files', metavar='sequence_files', nargs='+')
    parser.add_argument('--out-file', dest='destination_file',
            type=argparse.FileType('w'), default=sys.stdout,
            metavar='destination_file',
            help='Output destination. Default: STDOUT')
    parser.add_argument('--format', dest='output_format', default='tab',
        choices=('tab', 'csv', 'align'), help='''Specify output format
        as tab-delimited, CSV or aligned in a borderless table.
        Default is tab-delimited''')
    parser.add_argument('--width', dest='width', type=int, default=30,
         help='''Specify width of columns when --format is set to align.
         Defaults to a width of %(default)s.  Any output exceeding the column
         width will be truncated.''')


def _print_file_info(row, output_format, handle, width):
    """
    Write out information that describes a sequence file.
    """
    if 'tab' in output_format:
        handle.write("\t".join(row) + "\n")
    elif 'csv' in output_format:
        writer = csv.writer(handle, delimiter=',', quotechar='"',
                            quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(row)
    elif 'align' in output_format:
        row = map(lambda s: s.ljust(width)[:width], row)
        handle.write("".join(row) + "\n")


_HEADERS = ('name', 'alignment', 'min_len', 'max_len', 'avg_len',
              'num_seqs')
_SeqFileInfo = collections.namedtuple('SeqFileInfo', _HEADERS)


def _summarize_file(source_file):
    """
    Summarizes a sequence file, returning a tuple containing the name,
    whether the file is an alignment, minimum sequence length, maximum
    sequence length, average length, number of sequences.
    """
    is_alignment = True
    avg_length = None
    min_length = sys.maxint
    max_length = 0
    sequence_count = 0
    file_type = fileformat.lookup_file_type(
            os.path.splitext(source_file)[1])

    # Get an iterator and analyze the data.
    for record in SeqIO.parse(source_file, file_type):
        sequence_count += 1
        sequence_length = len(record)
        if max_length != 0:
            # If even one sequence is not the same length as the others,
            # we don't consider this an alignment.
            if sequence_length != max_length:
                is_alignment = False

        # Work on determining the length of the longest sequence.
        if sequence_length > max_length:
            max_length = sequence_length
        if sequence_length < min_length:
            min_length = sequence_length

        # Average length
        if sequence_count == 1:
            avg_length = float(sequence_length)
        else:
            avg_length = avg_length + ((sequence_length - avg_length) /
                                       sequence_count)

    # Handle an empty file:
    if avg_length is None:
        min_length = max_length = avg_length = 0

    return _SeqFileInfo(source_file, str(is_alignment).upper(),
            str(min_length), str(max_length), '{0:.2f}'.format(avg_length),
            str(sequence_count))

def action(arguments):
    """
    Given one more more sequence files, determine if the file is an
    alignment, the maximum sequence length and the total number of
    sequences.  Provides different output formats including tab
    (tab-delimited), csv and align (aligned as if part of a borderless
    table).
    """
    handle = arguments.destination_file
    output_format = arguments.output_format
    width = arguments.width

    # Create and write out the header row.
    header = ('name', 'alignment', 'min_len', 'max_len', 'avg_len',
              'num_seqs')

    _print_file_info(header, output_format=output_format,
                     handle=handle, width=width)

    with handle:
        # Go through all source files passed in, one by one.
        for source_file in arguments.source_files:
            row = _summarize_file(source_file)

            _print_file_info(row=row, output_format=output_format,
                             handle=handle, width=width)
