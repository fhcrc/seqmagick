"""
Common functions for subcommands
"""
import argparse
import copy
import functools
import os
import os.path

def cut_range(string):
    """
    A custom argparse 'type' to deal with sequences ranges such as 5:500.

    Returns a 0-based slice corresponding to the selection defined by the
    """
    value_range = string.split(':')
    if len(value_range) != 2:
        msg = "{0} is not a valid, 1-indexed range.".format(string)
        raise argparse.ArgumentTypeError(msg)

    start, stop = tuple(map(int, value_range))
    # Convert from 1-indexed to 0-indexed
    start -= 1
    if start < 0 or stop < start:
        msg = "{0} is not a valid, 1-indexed range.".format(string)
        raise argparse.ArgumentTypeError(msg)

    return slice(start, stop)


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


def partial_append_action(fn, argument_keys=None):
    """
    Creates a new class extending argparse.Action, which appends a
    partially-applied function to dest.

    The optional argument_keys argument should either be None (no additional
    arguments to fn) or an iterable of function keys to partially apply.
    """
    if isinstance(argument_keys, basestring):
        argument_keys = [argument_keys]
    argument_keys = argument_keys or []

    class PartialAppendAction(argparse.Action):
        def __init__(self,
                     option_strings,
                     dest,
                     const=None,
                     default=None,
                     required=False,
                     help=None,
                     type=None,
                     metavar=None,
                     nargs=None,
                     **kwargs):
            super(PartialAppendAction, self).__init__(
                option_strings=option_strings,
                dest=dest,
                nargs=len(argument_keys),
                const=const,
                default=default,
                required=required,
                metavar=metavar,
                type=type,
                help=help, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            items = copy.copy(getattr(namespace, self.dest, None)) or []

            # If no value was set default to empty list
            if values is None:
                values = []
            elif not isinstance(values, list):
                values = [values]

            if len(argument_keys) != len(values):
                raise ValueError("Unexpected number of values")

            # Generate keyword arguments for the input function
            kwargs = dict(zip(argument_keys, values))
            f = functools.partial(fn, **kwargs)
            functools.update_wrapper(f, fn)
            f.func_name = 'wrapped_' + fn.func_name
            items.append(f)
            setattr(namespace, self.dest, items)

    return PartialAppendAction


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
