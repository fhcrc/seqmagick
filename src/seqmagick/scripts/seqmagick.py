#! /usr/bin/env python

import os
import sys
import argparse

from cStringIO import StringIO
# Insert one level above project directory to path for testing.
sys.path.insert(0, "../..")
from seqmagick.magickwrap import MagickWrap
from seqmagick.fileformat import FileFormat



def main():
    action, arguments = parse_arguments()

    # We have nothing to do.
    if arguments is None:
        sys.exit(0)
  
    # Convert takes a single source, mogrify can accept multiple 
    # source files.
    sources = []
    destination = None
    if action == 'convert':
        sources = arguments.source_file
        destination = arguments.destination_file[0]
    elif action == 'mogrify':
        sources = arguments.source_files
        
    # Will want to catch control-C and exceptions here to remove 
    # destination file or temp file if script did not fully execute.
    if arguments is not None and action:
        wrap = MagickWrap(in_files=sources, out_file=destination, tmp_dir=arguments.tmp_dir)
        if action == 'muscle':
            wrap.create_muscle_alignment()
        if action == 'mogrify' or action == 'convert':
            wrap.transform(complement=arguments.complement,
                           cut=arguments.cut,
                           dashgap=arguments.dashgap, 
                           deduplicate_taxa=arguments.deduplicatetaxa,
                           deduplicate_sequences=arguments.deduplicateseqs,
                           degap=arguments.degap,
                           first_name_capture=arguments.firstname,
                           grep=arguments.grep,
                           lower=arguments.lower,
                           reverse=arguments.reverse,
                           upper=arguments.upper,
                          )

def parse_arguments():
    """
    Pull out command-line options and do different things, depending on the action specified.
    """

    argv = sys.argv[1:]

    # List of valid actions.
    actions = ('convert', 'mogrify', 'muscle', 'check', 'grep', 'head', 'sort', 'tail', 'help')

    # Create an argparse instance.
    parser = argparse.ArgumentParser(description='SeqMagick - Manipulate sequence files.')

    # Example command-line usage:
    #
    # seqmagick.py convert x.fasta y.phy
    # seqmagick.py mogrify --degap --upper --reverse x.fasta c.fasta t.fasta *.sth
    # seqmagick.py check x.fasta
    # seqmagick.py muscle x.fasta
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

    # Add arguments specific to the muscle action.
    if action in ('muscle'):
        parser.add_argument('source_file', type=sequence_file, nargs=1)
        parser.add_argument('destination_file', nargs=1)

    # Add arguments specific to the check action.
    if action in ('check'):
        parser.add_argument('--alphabet', dest='alphabet', help='To be implemented')
        pass


    # Add arguments shared between convert and mogrify
    if action in ('convert') or action in ('mogrify'):
        parser.add_argument('--cut', dest='cut', metavar="start:end", type=cut_range, 
                            help='Start and end positions for cutting sequences, : separated.  Includes last item.')
        parser.add_argument('--complement', action='store_true', help='Convert sequences into complements')
        parser.add_argument('--dashgap', action='store_true', help='Change . and : into - for all sequences')
        parser.add_argument('--deduplicateseqs', action='store_true', help='Remove any duplicate sequences by sequence content, keep the first instance seen')
        parser.add_argument('--deduplicatetaxa', action='store_true', help='Remove any duplicate sequences by ID, keep the first instance seen')
        parser.add_argument('--degap', action='store_true', help='Remove gaps in the sequence alignment')
        parser.add_argument('--firstname', action='store_true', help='Take only the first whitespace-delimited word as the name of the sequence') 
        parser.add_argument('--grep', dest='grep', help='Filter the sequences by regular expression in name') 
        parser.add_argument('--lower', action='store_true', help='Translate the sequences to lower case')
        parser.add_argument('--reverse', action='store_true', help='Reverse the order of sites in sequences')
        parser.add_argument('--strict', dest='data_type', metavar='data_type', 
                            help='Verify only IUPAC characters for "aa" or "nuc" are used')
        parser.add_argument('--translate', dest='destination_type', metavar='destination_type', 
                            help='Translate between amino acids and nucleotides, use "aa" or "nuc" as destination type')
        parser.add_argument('--upper', action='store_true', help='Translate the sequences to upper case')
        parser.add_argument('--wrap', action='store_true', help='')


    # Add arguments specific to the convert action.
    if action in ('convert'):
        parser.add_argument('source_file', type=sequence_file, nargs=1)
        parser.add_argument('destination_file', nargs=1, default=False)

    # Add arguments specific to the head action.
    if action in ('head'):
        #parser.add_argument('--', dest='', type=, help='')
        pass

    # Add arguments specific to the help action.
    if action in ('help'):
        print parser.description
        print_help(parser=parser)
        # Do not do anything in addition to printing help text, just exit.
        return action, None

    # Add arguments specific to the mogrify action.
    if action in ('mogrify'):
        parser.add_argument('source_files', type=sequence_file, nargs='+')

    # Add arguments specific to the sort action.
    if action in ('sort'):
        parser.add_argument('source_file', type=sequence_file, nargs=1)
        parser.add_argument('destination_file', nargs=1)

    # Add arguments specific to the tail action.
    if action in ('tail'):
        #parser.add_argument('--', dest='', type=, help='')
        pass

    # Add arguments common to all actions.
    if action in actions:
        parser.add_argument('--tmp', dest='tmp_dir', default='/tmp', help='Temporary directory for working file. Default is /tmp.')
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
    Print out general help and help for specific actions.  Does some string 
    manipulation to print text specific to an action, while continuing to use 
    argparse functionality.
    """
    if action is None:
        print """
SeqMagick actions include:
    check            Check integrity of a file.
    convert          Convert between sequence file formats and optionally,  
                     perform other operations.
    head             Print the top N records to STDOUT.
    mogrify          Perform in-place operations on a file containing sequences.
                     Can accept multiple source files.
    muscle           Create an alignment using muscle.
    sort             Sort sequences by length ascending, length descending,
                     name ascending or name descending.  Sorting is done 
                     in-memory.
    tail             Print the bottom N records to STDOUT.

See 'seqmagick.py help ACTION' for more information on a specific action.

"""

    if parser is not None and action:
        # Ugly, but it seems there is no clean way to put 'positional' arguments 
        # at the beginning of the usage section.
        help_text = parser.format_help()
        # Get rid of the action placeholder positional argument.
        help_text = help_text.replace('action source_file', 'source_file')
        help_text = help_text.replace('action [source_files [source_files ...]]', '[source_files [source_files ...]]')
        help_text = help_text.replace('  action\n', '')
        # Remove --help when an action is specified.
        help_text = help_text.replace('[-h]', action + '')
        # Have seen different amount of whitespace, will make sense to replace this with a regex.
        help_text = help_text.replace('  -h, --help        show this help message and exit\n', '')
        help_text = help_text.replace('  -h, --help            show this help message and exit\n', '')
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


