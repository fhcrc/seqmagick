#! /usr/bin/env python

import os
import sys
import argparse

# Insert one level above project directory to path for testing.
#sys.path.insert(0, "../..")
from seqmagick.magickwrap import MagickWrap
from seqmagick.fileformat import FileFormat


def main():
    action, arguments = parse_arguments()

    verbose = arguments.verbose
    debug = arguments.debug

    # Convert takes a single source, mogrify can accept multiple
    # source files.
    sources = []
    destination = None

    if action == 'convert' or action == 'muscle':
        sources = arguments.source_file
        destination = arguments.destination_file[0]
    elif action == 'info':
        sources = arguments.source_files
        destination = arguments.destination_file
    elif action == 'mogrify':
        sources = arguments.source_files

    if verbose: print 'Arguments parsed.  Action to be performed: ' + action
    if debug: print 'DEBUG: command-line arguments:\n' + str(arguments)

    # Will want to catch control-C and exceptions here to remove
    # destination file or temp file if script did not fully execute.
    if arguments is not None and action:
        if verbose: print 'Creating MagickWrap instance.'
        wrap = MagickWrap(in_files=sources, out_file=destination,
                          tmp_dir=arguments.tmp_dir, debug=debug,
                          verbose=verbose)
        if action == 'info':
            if verbose: print 'Performing info action.'
            wrap.describe_sequence_files(output_format=arguments.format,
                                         width=arguments.width,
                                        )
        if action == 'muscle':
            if verbose: print 'Performing muscle alignment.'
            wrap.create_muscle_alignment()
        if action == 'mogrify' or action == 'convert':
            if verbose: print 'Performing ' + action + '.'
            wrap.transform(cut=arguments.cut,
                           dashgap=arguments.dashgap,
                           deduplicate_taxa=arguments.deduplicatetaxa,
                           deduplicate_sequences=arguments.deduplicateseqs,
                           ungap=arguments.ungap,
                           first_name_capture=arguments.firstname,
                           head=arguments.head,
                           input_format=arguments.input_format,
                           linewrap=arguments.linewrap,
                           lower=arguments.lower,
                           max_length=arguments.max_length,
                           min_length=arguments.min_length,
                           name_suffix=arguments.name_suffix,
                           name_prefix=arguments.name_prefix,
                           output_format=arguments.output_format,
                           pattern_include=arguments.pattern_include,
                           pattern_exclude=arguments.pattern_exclude,
                           reverse=arguments.reverse,
                           reverse_complement=arguments.reverse_complement,
                           sort=arguments.sort,
                           strip_range=arguments.strip_range,
                           squeeze=arguments.squeeze,
                           tail=arguments.tail,
                           transcribe=arguments.transcribe,
                           translate=arguments.translate,
                           upper=arguments.upper,
                           )


def parse_arguments(action_arguments=None):
    """
    Extract command-line arguments for different actions.
    """
    parser = argparse.ArgumentParser(description='SeqMagick - Manipulate sequence files.', prog='seqmagick')

    # Setup sub-commands
    subparsers = parser.add_subparsers(dest='subparser_name')
    # Check
    parser_check = subparsers.add_parser('check', help='Check integrity of a file.')
    # To be implemented
    #parser_check.add_argument('--alphabet', dest='alphabet', help='To be implemented')
    parser_check.add_argument('source_file', type=sequence_file, nargs=1)
    # Convert
    parser_convert = add_arguments(subparsers.add_parser('convert',
                                   help='Convert between sequence file formats and optionally, '
                                   'perform other operations'))
    parser_convert.add_argument('source_file', type=sequence_file, nargs=1)
    parser_convert.add_argument('destination_file', nargs=1, default=False)
    # Help
    parser_help = subparsers.add_parser('help', help='Help for actions')
    parser_help.add_argument('action', nargs=1)
    # Info
    parser_info = subparsers.add_parser('info', help='Describe sequence files.')
    parser_info.add_argument('source_files', metavar='sequence_files',
                             type=sequence_file, nargs='+')
    parser_info.add_argument('--outfile', dest='destination_file', metavar='destination_file',
                             help='Write output to a file.  If unspecified, ' +
                             'defaults to STDOUT')
    parser_info.add_argument('--format', dest='format', default='tab',
                             choices=['tab', 'csv', 'align'],
                             help='Specify output format as tab-delimited, CSV or ' +
                             'aligned in a borderless table.  Default is tab-delimited')
    parser_info.add_argument('--width', dest='width', type=int, default=30,
                             help='Specify width of columns when --format is set to align. ' +
                             'Defaults to a width of 30.  Any output exceeding the column width ' +
                             'will be truncated ')
    # Mogrify
    parser_mogrify = add_arguments(subparsers.add_parser('mogrify',
                                   help='Perform in-place operations on a file '
                                   ' containing sequences.  Can accept multiple source files. '))
    parser_mogrify.add_argument('source_files', type=sequence_file, nargs='+')
    # Muscle
    parser_muscle = subparsers.add_parser('muscle', help=' Create an alignment using muscle.')
    parser_muscle.add_argument('source_file', type=sequence_file, nargs=1)
    parser_muscle.add_argument('destination_file', nargs=1)

    # With the exception of 'help', all subcommands share a certain
    # number of arguments, which are added here.
    for subcommand in subparsers.choices.keys():
        if subcommand == 'help': continue
        subparser = subparsers.choices[subcommand]
        subparser.add_argument('--tmp', dest='tmp_dir', default='/tmp',
                                help='Temporary directory for working file. Default is /tmp.')
        subparser.add_argument('--debug', action='store_true', help='Enable debug output')
        subparser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    # Override arguments passed at the command-line to support seqmagick help <action>
    if action_arguments:
        arguments = parser.parse_args(action_arguments)
    else:
        arguments = parser.parse_args()

    action = arguments.subparser_name

    # Support seqmagick help <action> by simply having this function call itself and
    # translate the arguments into something that argparse can work with.
    if (action == 'help'):
        parse_arguments(action_arguments=[str(arguments.action[0]), '-h'])

    return action, arguments


def add_arguments(subparser):
    """
    Add shared arguments to the convert or mogrify subparser.
    """
    subparser.add_argument('--cut', dest='cut', metavar="start:end", type=cut_range,
                        help='1-indexed start and end positions for cutting sequences, : separated.  Includes last item.')
    subparser.add_argument('--dashgap', action='store_true',
                        help='Change . and : into - for all sequences')
    subparser.add_argument('--deduplicateseqs', action='store_true',
                        help='Remove any duplicate sequences by sequence content, keep the first instance seen')
    subparser.add_argument('--deduplicatetaxa', action='store_true',
                        help='Remove any duplicate sequences by ID, keep the first instance seen')
    subparser.add_argument('--firstname', action='store_true',
                        help='Take only the first whitespace-delimited word as the name of the sequence')
    subparser.add_argument('--head', metavar='N', dest='head', type=int, help='Trim down to top N sequences')
    subparser.add_argument('--linewrap', dest='linewrap', metavar='N', type=int,
            help='Adjust line wrap for sequence strings.  When N is 0, '
                 'all line breaks are removed. Only fasta files are supported for the output format.')
    subparser.add_argument('--lower', action='store_true', help='Translate the sequences to lower case')
    subparser.add_argument('--maxlength', dest='max_length', metavar='N',
                           type=int, help='Discard any sequences beyond '
                           'the specified maximum length.  This operation occurs *before* all length-changing '
                           'options such as cut and squeeze.')
    subparser.add_argument('--minlength', dest='min_length', metavar='N',
                           type=int, help='Discard any sequences less than '
                           'the specified minimum length.  This operation occurs *before* all length-changing '
                           'options such as cut and squeeze.')
    subparser.add_argument('--name-suffix', metavar='SUFFIX', dest='name_suffix',
                        help='Append a suffix to all IDs.')
    subparser.add_argument('--name-prefix', metavar='PREFIX', dest='name_prefix',
                        help='Insert a prefix for all IDs.')
    subparser.add_argument('--pattern-include', metavar='regex', dest='pattern_include',
                        help='Filter the sequences by regular expression in name')
    subparser.add_argument('--pattern-exclude', metavar='regex', dest='pattern_exclude',
                        help='Filter out sequences by regular expression in name')
    subparser.add_argument('--reverse', action='store_true', help='Reverse the order of sites in sequences')
    subparser.add_argument('--reversecomplement', dest='reverse_complement', action='store_true',
                        help='Convert sequences into reverse complements')
    subparser.add_argument('--sort', dest='sort',
                        choices=['length-asc', 'length-desc', 'name-asc', 'name-desc'],
                        help='Perform sorting by length or name, ascending or descending.  '
                        'ASCII sorting is performed for names')
    subparser.add_argument('--striprange', dest='strip_range', action='store_true',
                        help='Strip ranges from sequences IDs, matching </x-y>')

    subparser.add_argument('--squeeze', action='store_true',
                        help='Remove any gaps that are present in the same position across ' +
                        'all sequences in an alignment')
    subparser.add_argument('--tail', metavar='N', dest='tail', type=int, help='Trim down to bottom N sequences')
    # To be implemented
    #subparser.add_argument('--translate', dest='destination_type', metavar='destination_type',
    #                    help='Translate between amino acids and nucleotides, use "aa" or "nuc" as destination type')
    subparser.add_argument('--transcribe', dest='transcribe',
                        choices=['dna2rna', 'rna2dna'],
                        help='Transcription and back transcription for generic DNA and RNA.  '
                        'Source sequences must be the correct alphabet or this action will likely '
                        'produce incorrect results.')
    subparser.add_argument('--translate', dest='translate',
                        choices=['dna2protein', 'rna2protein', 'dna2proteinstop', 'rna2proteinstop'],
                        help='Translate from generic DNA/RNA to proteins. Options with '
                        '"stop" suffix will NOT translate through stop codons .'
                        'Source sequences must be the correct alphabet or this action will likely '
                        'produce incorrect results.')
    subparser.add_argument('--ungap', action='store_true', help='Remove gaps in the sequence alignment')
    subparser.add_argument('--upper', action='store_true', help='Translate the sequences to upper case')

    format_group = subparser.add_argument_group('Format Options')
    format_group.add_argument('--input-format', metavar='Format',
            help="Input file format (default: determine from extension")
    format_group.add_argument('--output-format', metavar='Format',
            help="Output file format (default: determine from extension)")

    return subparser


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


def sequence_file(sequence_file):
    """
    A custom argparse 'type' to make sure sequence files exist and the type is supported.
    """
    if os.access(sequence_file, os.R_OK) and os.path.isfile(sequence_file):
        return sequence_file
    else:
        raise IOError(sequence_file + " not found or is not readable.")


if __name__ == '__main__':
    sys.exit(main())

