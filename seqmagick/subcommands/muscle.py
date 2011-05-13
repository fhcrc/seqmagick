"""
Muscle Action
"""

import logging
import shlex
import subprocess

from Bio.Align.Applications import MuscleCommandline

import common


def build_parser(parser):
    parser.add_argument('source_file', type=common.sequence_file)
    parser.add_argument('destination_file')


def action(arguments):
    """
    Use BioPython muscle wrapper to create an alignment.
    """
    muscle_command = MuscleCommandline(
            input=arguments.source_file, out=arguments.destination_file)
    logging.debug('Muscle command: %s', muscle_command)

    return_code = subprocess.check_call(shlex.split(str(muscle_command)))
    return return_code
