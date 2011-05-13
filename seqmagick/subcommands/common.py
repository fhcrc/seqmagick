"""
Common functions for subcommands
"""
import argparse
import os
import os.path

def cut_range(string):
    """
    A custom argparse 'type' to deal with sequences ranges such as 5:500.
    """
    value_range = string.split(':')
    # We want integers, and something immutable.
    value_range = tuple(map(int, value_range))

    # Make sure the value range looks sane.
    if (len(value_range) != 2 or value_range[0] < 1 or
            value_range[0] > value_range[1]):
        msg = "%s is not a valid, 1-indexed range." % string
        raise argparse.ArgumentTypeError(msg)
    return value_range


def sequence_file(sequence_file):
    """
    A custom argparse 'type' to make sure sequence files exist and the type is
    supported.
    """
    if os.access(sequence_file, os.R_OK) and os.path.isfile(sequence_file):
        return sequence_file
    else:
        raise argparse.ArgumentTypeError(sequence_file +
                " not found or is not readable.")
