"""
Common functions for subcommands
"""
import argparse
import functools
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
        msg = "{0} is not a valid, 1-indexed range.".format(string)
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

def typed_range(type_func, minimum, maximum):
    """
    Require variables to be of the specified type, between minimum and maximum
    """
    @functools.wraps(type_func)
    def inner(string):
        result = type_func(string)
        if not result >= minimum and result <= maximum:
            raise argparse.ArgumentTypeError(
                    "Please provide a value between {0} and {1}".format(
                        minimum, maximum))
        return result
    return inner


def positive_value(target_type):
    """
    Wraps target_type in a function that requires the parsed argument
    be >= 0
    """
    def inner(string):
        value = target_type(string)
        if not value >= 0:
            raise argparse.ArgumentTypeError("Invalid positive number: " +
                    string)
        return value

    return inner
