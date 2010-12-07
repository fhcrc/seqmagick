#! /usr/bin/env python

import os
import sys
import argparse

from cStringIO import StringIO
from SeqMagick.magickwrap import MagickWrap
from SeqMagick.fileformat import FileFormat

# Append top-level project directory to path for testing.
#sys.path.append(".")
#sys.path.append("..")
#from magickwrap import MagickWrap 
#from fileformat import FileFormat






def main():
    action, arguments = parse_arguments()

    if arguments is not None and action is not None:
        wrap = MagickWrap(in_file=arguments.source_file[0], out_file=arguments.destination_file[0])
        if action == 'convert':
            wrap.convert_format()

def parse_arguments():
    """
    Pull out command-line options and do different things, depending on the action specified.
    """

    argv = sys.argv[1:]

    # List of valid actions.
    actions = ('convert', 'mogrify', 'align', 'check', 'grep', 'head', 'tail', 'help')

    # Create an argparse instance.
    parser = argparse.ArgumentParser(description='SeqMagick - Manipulate sequence files.')

    # Example command-line usage:
    #
    # seqmagick.py convert x.fasta y.phy
    # seqmagick.py mogrify --degap --upper --reverse x.fasta 
    # seqmagick.py check x.fasta
    # seqmagick.py align x.fasta
    # seqmagick.py head -20 x.fasta
    # seqmagick.py head -20 x.fasta
    # seqmagick.py tail -20 x.fasta
    # seqmagick.py help mogrify

    # An action is required, print help if nothing specified.
    if len(argv) == 0:
        print parser.description
        print_help(parser=parser)
        return None, None

    # Determine the action.
    action = 'help' # Default to help action.
    if len(argv) > 1:
        # Try to grab help actions first (second argument), 
        # then attempt to grab actions from first argument. 
        if argv[1].lower() in actions:
            action = argv[1].lower()
        elif argv[0].lower() in actions:
             action = argv[0].lower()

    # Grab first argument which will usually be the action.
    # If positional...
    parser.add_argument('action')

    # Check to see if help should be printed.  Do not execute 
    # any commands when help is to be displayed, make sure 
    # action is set to 'help'
    if '-h' in argv or '--help' in argv:
        action = 'help' 

    # Argument groups might be an alternative to the way things are done below.

    # Add arguments specific to the align action.
    if action in ('align'):
        #parser.add_argument('--', dest='', type=, help='')
        pass

    # Add arguments specific to the align action.
    if action in ('check'):
        #parser.add_argument('--', dest='', type=, help='')
        pass

    # Add arguments specific to the align action.
    if action in ('convert'):
        parser.add_argument('source_file', type=sequence_file, nargs=1)
        parser.add_argument('destination_file', nargs=1)

    # Add arguments specific to the grep action.
    if action in ('grep'):
        #parser.add_argument('--', dest='', type=, help='')
        pass

    # Add arguments specific to the head action.
    if action in ('head'):
        #parser.add_argument('--', dest='', type=, help='')
        pass

    # Add arguments specific to the help action.
    if action in ('help'):
        print_help(parser=parser)
        # Do not do anything in addition to printing help text, just exit.
        return action, None

    # Add arguments specific to the mogrify action.
    if action in ('mogrify'):
        parser.add_argument('--cut', dest='cut', metavar="start:end", type=cut_range, 
                            help='Start and end positions for cutting sequences, : separated')
        parser.add_argument('--dashgap', action='store_true', help='Change . and : into - for all sequences')
        parser.add_argument('--degap', action='store_true', help='Remove gaps in the sequence alignment')
        parser.add_argument('--lower', action='store_true', help='Translate the sequences to lower case')
        parser.add_argument('--reverse', dest='', help='')
        parser.add_argument('--sort', dest='', help='')
        parser.add_argument('--strict', dest='', help='')
        parser.add_argument('--translate', dest='', help='')
        parser.add_argument('--upper', dest='', help='Translate the sequences to upper case.')
        parser.add_argument('--wrap', dest='', help='')
        #parser.add_argument('--', dest='', type=, help='')

    # Add arguments specific to the tail action.
    if action in ('tail'):
        #parser.add_argument('--', dest='', type=, help='')
        pass

    # Add arguments common to all actions.
    if action in actions:
        parser.add_argument('--debug', dest='', help='')
        parser.add_argument('--verbose', dest='', help='')

    # Print help relevant to the action.
    if len(argv) == 2 and argv[0].lower() == 'help' and argv[1].lower() in actions:
        print_help(action=argv[1].lower(), parser=parser)
        return action, None

    #  Try to return action and arguments.  Capture STDERR and report any 
    #  errors that are detected.
    capture = StringIO()
    stderr_orig = sys.stderr
    try:
        sys.stderr = capture
        arguments = parser.parse_args()
        sys.stderr = stderr_orig
        return action, arguments
    except:
        sys.stderr = stderr_orig
        print_help(action=action, parser=parser, message=capture.getvalue())
        raise
    finally:
        capture.close()

def print_help(action=None, parser=None, message=None):
    """
    """
    if action is None:
        print """
SeqMagick actions include:
    align            Create an alignment.
    check            Check integrity of a file.
    convert          Convert between sequence file formats and optionally,  
                     perform other operations.
    grep             Filter the sequences by regular expression in name 
                     to STDOUT.
    head             Print the top N records to STDOUT.
    mogrify          Perform in-place operations on a file containing sequences.
    tail             Print the bottom N records to STDOUT.

See 'seqmagick.py help ACTION' for more information on a specific action.

"""

    if parser is not None and action is not None:
        # Ugly, but it seems there is no clean way to put 'positional' arguments 
        # at the beginning of the usage section.
        help_text = parser.format_help()
        # Get rid of the action placeholder positional argument.
        help_text = help_text.replace('action source_file', 'source_file')
        help_text = help_text.replace('  action\n', '')
        # Remove --help when an action is specified.
        help_text = help_text.replace('[-h]', action + '')
        help_text = help_text.replace('  -h, --help        show this help message and exit\n', '')
        help_text = help_text.replace('optional arguments', 'optional arguments for the ' + action + ' action')
        print help_text

    # Print just the error message if one came through, which should be on 
    # the second to last line.
    if message is not None:
        lines = message.split("\n")
        last_line = lines[len(lines) - 2]
        if ('error' in last_line.lower()):
            print last_line

# Custom argparse types.

def sequence_file(sequence_file):
    """
    A custom argparse 'type' to make sure sequence files exist and the type is supported.
    """
    FileFormat.lookup_file_type(os.path.splitext(sequence_file)[1])
    if os.access(sequence_file, os.R_OK) and os.path.isfile(sequence_file):
        return sequence_file
    else:
        raise Exception, sequence_file + " not found or is not readable."
    

def cut_range(string):
    """
    A custom argparse 'type' to deal with sequences ranges such as 5:500.
    """
    value_range = string.split(':')
    # We want integers, and something immutable.
    value_range = tuple(map(int, value_range))

    # Make sure the value range looks sane.
    if len(value_range) != 2 or value_range[0] < 1 or  value_range[0] >  value_range[1]:
        msg = "%s is not a valid, 1-indexed range." % string
        raise argparse.ArgumentTypeError(msg)
    return value_range


if __name__ == '__main__':
    sys.exit(main())


