"""
Modify sequence file(s) in place.
"""

import convert
import logging
import os.path
import shutil
import tempfile


def build_parser(parser):
    """
    """
    convert.add_options(parser)

    parser.add_argument('input_files', metavar="sequence_file", nargs='+',
            help="Sequence file(s) to mogrify")

    return parser


def action(arguments):
    """
    Run mogrify.  Most of the action is in convert, this just creates a temp
    file for the output.
    """
    for input_file in arguments.input_files:
        logging.info(input_file)
        bn = os.path.basename(input_file)
        # Generate a temporary file
        with tempfile.NamedTemporaryFile(prefix='smagick', suffix=bn,
                delete=False) as tf:
            temp_name = tf.name

        convert.transform_file(input_file, temp_name, arguments)
        # Overwrite the original file
        shutil.move(temp_name, input_file)

