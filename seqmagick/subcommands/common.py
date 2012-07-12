"""
Common functions for subcommands
"""
import argparse
import contextlib
import copy
import functools
import os
import os.path
import signal
import sys
import tempfile

def get_umask():
    """
    Gets the current umask
    """
    current_umask = os.umask(0777)
    os.umask(current_umask)
    return current_umask

def apply_umask(permission=0666, umask=None):
    """
    Masks the provided permission with a umask.

    If umask is not given, the current umask is used.
    """
    if umask is None:
        umask = get_umask()
    return permission & (~umask)

@contextlib.contextmanager
def atomic_write(path, permissions=None, **kwargs):
    """
    Open a file for atomic writing.

    Generates a temp file, renames to value of ``path``.

    Additional arguments are passed to tempfile.NamedTemporaryFile
    """
    if permissions is None:
        permissions = apply_umask()
    # Handle stdout:
    if path == '-':
        yield sys.stdout
    else:
        base_dir = os.path.dirname(path)
        kwargs['suffix'] = '.' + os.path.splitext(path)[1]
        tf = tempfile.NamedTemporaryFile(dir=base_dir, delete=False,
                                         **kwargs)
        try:
            with tf:
                yield tf
            # Move
            os.rename(tf.name, path)
            os.chmod(path, permissions)
        except:
            os.remove(tf.name)
            raise

def sequence_slices(string):
    """
    Parses a list of slices from a string of format:

    start1:end1[,start2:end2[,start2:end3]] etc
    """
    slices = string.split(',')
    return [cut_range(i) for i in slices]

def cut_range(string):
    """
    A custom argparse 'type' to deal with sequences ranges such as 5:500.

    Returns a 0-based slice corresponding to the selection defined by the slice
    """
    value_range = string.split(':')
    if len(value_range) == 1:
        start = int(value_range[0])
        stop = start
    elif len(value_range) == 2:
        start, stop = tuple(int(i) if i else None for i in value_range)
    else:
        msg = "{0} is not a valid, 1-indexed range.".format(string)
        raise argparse.ArgumentTypeError(msg)

    # Convert from 1-indexed to 0-indexed
    if start is not None:
        start -= 1

    if ((start or 0) < 0
            or (stop if stop is not None else sys.maxint) < (start or 0)):
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

def _exit_on_signal(sig, status=None, message=None):
    def exit(sig, frame):
        if message:
            print >> sys.stderr, message
        raise SystemExit(status)
    signal.signal(sig, exit)

def exit_on_sigint(status=1, message="Canceled."):
    """
    Set program to exit on SIGINT, with provided status and message.
    """
    _exit_on_signal(signal.SIGINT, status, message)

def exit_on_sigpipe(status=None):
    """
    Set program to exit on SIGPIPE
    """
    _exit_on_signal(signal.SIGPIPE, status)
