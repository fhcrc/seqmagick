"""
Modify sequence file(s) in place.
"""

import argparse
import logging
import os.path

from . import convert, common

def build_parser(parser):
    """
    """
    convert.add_options(parser)

    parser.add_argument(
        'input_files', metavar="sequence_file", nargs='+',
        type=argparse.FileType('r+'),
        help="Sequence file(s) to mogrify")

    return parser


def action(arguments):
    """
    Run mogrify.  Most of the action is in convert, this just creates a temp
    file for the output.
    """
    for input_file in arguments.input_files:
        logging.info(input_file)
        # Generate a temporary file
        with common.atomic_write(input_file.name) as tf:
            convert.transform_file(input_file, tf, arguments)
