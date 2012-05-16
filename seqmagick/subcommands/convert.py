"""
Convert between sequence formats
"""
import argparse
import functools
import logging
import os.path

from Bio import Alphabet, SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqIO import FastaIO
from seqmagick import transform
from seqmagick.fileformat import from_extension

from . import common

ALPHABETS = {
        'dna': Alphabet.generic_dna,
        'dna-ambiguous': IUPAC.ambiguous_dna,
        'protein': Alphabet.generic_protein,
        'rna': Alphabet.generic_rna,
        'rna-ambiguous': IUPAC.ambiguous_rna,
}

def add_options(parser):
    """
    Add optional arguments to the parser
    """
    partial_action = common.partial_append_action
    file_mods = parser.add_argument_group("Sequence File Modification")
    file_mods.add_argument('--line-wrap', dest='line_wrap', metavar='N',
        type=int, help='Adjust line wrap for sequence strings.  '
        'When N is 0, all line breaks are removed. Only fasta files '
        'are supported for the output format.')
    file_mods.add_argument('--sort', dest='sort',
        choices=['length-asc', 'length-desc', 'name-asc', 'name-desc'],
        help='Perform sorting by length or name, ascending or descending. '
        'ASCII sorting is performed for names')

    seq_mods = parser.add_argument_group("Sequence Modificaton")
    seq_mods.add_argument('--apply-function', type=module_function,
            metavar='/path/to/module.py:function_name[:parameter]',
            help="""Specify a custom function to apply to the input sequences,
            specified as /path/to/file.py:function_name. Function should accept
            an iterable of Bio.SeqRecord objects, and yield SeqRecords. If the
            parameter is specified, it will be passed as a string as the second
            argument to the function. Specify more than one to chain.""",
            default=[], action='append')
    seq_mods.add_argument('--cut', dest='transforms',
            metavar="start:end[,start2:end2]",
            type=common.sequence_slices,
            action=partial_action(transform.multi_cut_sequences, 'slices'),
            help="""Keep only the residues within the 1-indexed start and end
            positions specified, : separated. Includes last item. Start or end
            can be left unspecified to indicate start/end of sequence.""")
    seq_mods.add_argument('--relative-to', dest='cut_relative', metavar='ID',
            help="""Apply --cut relative to the indexes of non-gap residues in
            sequence identified by ID""")
    seq_mods.add_argument('--drop', dest='transforms',
            metavar='start:end[,start2:end2]',
            type=common.sequence_slices,
            action=partial_action(transform.drop_columns, 'slices'),
            help="""Remove the residues at the specified indices. Same format as `--cut`.""")
    seq_mods.add_argument('--dash-gap',
            action=partial_action(transform.dashes_cleanup), dest='transforms',
            help="""Replace any of the characters "?.:~" with a "-" for all
            sequences""")
    seq_mods.add_argument('--lower',
            action=partial_action(transform.lower_sequences),
            dest='transforms', help='Translate the sequences to lower case')
    seq_mods.add_argument('--mask', metavar="start1:end1[,start2:end2]",
            action=partial_action(transform.multi_mask_sequences, 'slices'),
            type=common.sequence_slices, dest='transforms', help="""Replace
            residues in 1-indexed slice with gap-characters. If --relative-to
            is also specified, coordinates are relative to the sequence ID
            provided.""")
    seq_mods.add_argument('--reverse',
            action=partial_action(transform.reverse_sequences),
            dest='transforms', help='Reverse the order of sites in sequences')
    seq_mods.add_argument('--reverse-complement', dest='transforms',
            action=partial_action(transform.reverse_complement_sequences),
            help='Convert sequences into reverse complements')
    seq_mods.add_argument('--squeeze', action=partial_action(transform.squeeze),
            dest='transforms',
            help='''Remove any gaps that are present in the same
            position across all sequences in an alignment (equivalent to
            --squeeze-threshold=1.0)''')
    seq_mods.add_argument('--squeeze-threshold', dest='transforms',
            action=partial_action(transform.squeeze, 'gap_threshold'),
            type=common.typed_range(float, 0.0, 1.0),
            metavar='PROP', help="""Trim columns from an alignment which
            have gaps in least the specified proportion of sequences.""")
    seq_mods.add_argument('--transcribe', dest='transforms',
            action=partial_action(transform.transcribe, 'transcribe'),
            choices=('dna2rna', 'rna2dna'), help="""Transcription and back
            transcription for generic DNA and RNA. Source sequences must be the
            correct alphabet or this action will likely produce incorrect
            results.""")
    seq_mods.add_argument('--translate', dest='transforms',
            action=partial_action(transform.translate, 'translate'),
            choices=['dna2protein', 'rna2protein', 'dna2proteinstop',
                'rna2proteinstop'], help="""Translate from generic DNA/RNA to
            proteins. Options with "stop" suffix will NOT translate through
            stop codons .  Source sequences must be the correct alphabet or
            this action will likely produce incorrect results.""")
    seq_mods.add_argument('--ungap',
            action=partial_action(transform.ungap_sequences),
            dest='transforms', help='Remove gaps in the sequence alignment')
    seq_mods.add_argument('--upper',
            action=partial_action(transform.upper_sequences),
            dest='transforms', help='Translate the sequences to upper case')

    seq_select = parser.add_argument_group("Record Selection")

    seq_select.add_argument('--deduplicate-sequences',
        action='store_const', const=None, default=False,
         dest='deduplicate_sequences', help='Remove any duplicate sequences '
         'by sequence content, keep the first instance seen')
    seq_select.add_argument('--deduplicated-sequences-file', action='store',
        metavar='FILE', dest='deduplicate_sequences', default=False,
        type=argparse.FileType('w'),
        help='Write all of the deduplicated sequences to a file')
    seq_select.add_argument('--deduplicate-taxa',
            action=partial_action(transform.deduplicate_taxa),
            dest='transforms', help="""Remove any duplicate sequences by ID,
            keep the first instance seen""")
    seq_select.add_argument('--exclude-from-file', metavar='FILE',
            type=argparse.FileType('r'), help="""Filter sequences, removing
            those sequence IDs in the specified file""", dest='transforms',
            action=partial_action(transform.exclude_from_file, 'handle'))
    seq_select.add_argument('--include-from-file', metavar='FILE',
            type=argparse.FileType('r'), help="""Filter sequences, keeping only
            those sequence IDs in the specified file""", dest='transforms',
            action=partial_action(transform.include_from_file, 'handle'))
    seq_select.add_argument('--head', metavar='N', dest='transforms', type=int,
            action=partial_action(transform.head, 'head'), help="""Trim
            down to top N sequences""")
    seq_select.add_argument('--max-length', dest='transforms', metavar='N',
            action=partial_action(transform.max_length_discard, 'max_length'),
            type=int, help="""Discard any sequences beyond the specified
            maximum length.  This operation occurs *before* all length-changing
            options such as cut and squeeze.""")
    seq_select.add_argument('--min-length', dest='transforms', metavar='N',
            action=partial_action(transform.min_length_discard, 'min_length'),
            type=int, help="""Discard any sequences less than the specified
            minimum length.  This operation occurs *before* cut and squeeze.""")
    seq_select.add_argument('--min-ungapped-length', metavar='N',
            action=partial_action(transform.min_ungap_length_discard,
                'min_length'), type=int, help="""Discard any sequences less
                than the specified minimum length, excluding gaps. This
                operation occurs *before* cut and squeeze.""",
                dest='transforms')
    seq_select.add_argument('--pattern-include', metavar='regex',
            action=partial_action(transform.name_include, 'filter_regex'),
            dest='transforms', help="""Filter the sequences by regular
            expression in name""")
    seq_select.add_argument('--pattern-exclude', metavar='regex',
            action=partial_action(transform.name_exclude, 'filter_regex'),
            dest='transforms', help="""Filter the sequences by regular
            expression in name""")
    seq_select.add_argument('--prune-empty',
            action=partial_action(transform.prune_empty), dest='transforms',
            help="Prune sequences containing only gaps ('-')")
    seq_select.add_argument('--seq-pattern-include', metavar='regex',
            action=partial_action(transform.seq_include, 'filter_regex'),
            dest='transforms', help="""Filter the sequences by regular
            expression in sequence""")
    seq_select.add_argument('--seq-pattern-exclude', metavar='regex',
            action=partial_action(transform.seq_exclude, 'filter_regex'),
            dest='transforms', help="""Filter the sequences by regular
            expression in sequence""")
    seq_select.add_argument('--tail', metavar='N', dest='transforms', type=int,
            action=partial_action(transform.tail, 'tail'),
        help='Trim down to bottom N sequences')

    id_mods = parser.add_argument_group("Sequence ID Modification")
    id_mods.add_argument('--first-name',
            action=partial_action(transform.first_name_capture),
            dest='transforms', help='''Take only the first whitespace-delimited
            word as the name of the sequence''')
    id_mods.add_argument('--name-suffix', metavar='SUFFIX',
            action=partial_action(transform.name_append_suffix, 'suffix'),
            dest='transforms', help='Append a suffix to all IDs.')
    id_mods.add_argument('--name-prefix', metavar='PREFIX',
            action=partial_action(transform.name_insert_prefix, 'prefix'),
            dest='transforms', help="""Insert a prefix for all
            IDs.""")
    id_mods.add_argument('--pattern-replace', nargs=2,
            metavar=('search_pattern', 'replace_pattern'),
            action=partial_action(transform.name_replace, ('search_regex',
                'replace_pattern')),
            dest='transforms', help="""Replace regex pattern "search_pattern"
            with "replace_pattern" in sequence ID""")
    id_mods.add_argument('--strip-range', dest='transforms',
            action=partial_action(transform.strip_range), help="""Strip ranges
            from sequences IDs, matching </x-y>""")

    format_group = parser.add_argument_group('Format Options')
    format_group.add_argument('--input-format', metavar='Format',
            help="Input file format (default: determine from extension)")
    format_group.add_argument('--output-format', metavar='Format',
            help="Output file format (default: determine from extension)")

    parser.add_argument('--alphabet', choices=ALPHABETS,
            help="""Input alphabet. Required for writing NEXUS.""")

    return parser


def build_parser(parser):
    """
    Add shared arguments to the convert or mogrify parser.
    """
    add_options(parser)
    parser.add_argument('source_file', type=argparse.FileType('r'),
                        help="Input sequence file")
    parser.add_argument('dest_file', help="Output file")

    return parser

def transform_file(source_file, destination_file, arguments):
    # Get just the file name, useful for naming the temporary file.
    file_ext = os.path.splitext(source_file.name)[1]
    source_file_type = (arguments.input_format or from_extension(file_ext))

    output_ext = os.path.splitext(getattr(destination_file, 'name', ''))[1]

    destination_file_type = (arguments.output_format or
            from_extension(output_ext))

    # Get an iterator.
    sorters = {'length': transform.sort_length,
               'name': transform.sort_name,}
    directions = {'asc': 1, 'desc': 0}
    if arguments.sort:
        # Sorted iterator
        key, direction = arguments.sort.split('-')
        records = sorters[key](source_file=source_file,
                source_file_type=source_file_type,
                direction=directions[direction])
    else:
        # Unsorted iterator.
        records = SeqIO.parse(source_file, source_file_type,
                alphabet=ALPHABETS.get(arguments.alphabet))


    #########################################
    # Apply generator functions to iterator.#
    #########################################

    # Apply all the transform functions in transforms
    if arguments.transforms:

        # Special case handling for --cut and --relative-to
        if arguments.cut_relative:
            for o, n in ((transform.multi_cut_sequences,
                          transform.cut_sequences_relative),
                         (transform.multi_mask_sequences,
                          transform.mask_sequences_relative)):
                # Add a function to trim any columns which are gaps in the
                # sequence ID
                try:
                    f = next(f for f in arguments.transforms
                             if f.func == o)
                except StopIteration:
                    continue
                i = arguments.transforms.index(f)
                arguments.transforms.pop(i)
                arguments.transforms.insert(i,
                        functools.partial(n,
                            record_id=arguments.cut_relative, **f.keywords))

        for function in arguments.transforms:
            records = function(records)

    if (arguments.deduplicate_sequences or
            arguments.deduplicate_sequences is None):
        records = transform.deduplicate_sequences(
            records, arguments.deduplicate_sequences)

    # Apply all the partial functions
    if arguments.apply_function:
        for apply_function in arguments.apply_function:
            records = apply_function(records)

    # Only the fasta format is supported, as SeqIO.write does not have a 'wrap'
    # parameter.
    if (arguments.line_wrap is not None and destination_file_type == 'fasta'
            and source_file_type == 'fasta'):
        logging.info("Attempting to write fasta with %d line breaks.",
                arguments.line_wrap)

        with destination_file:
            writer = FastaIO.FastaWriter(
                destination_file, wrap=arguments.line_wrap)
            writer.write_file(records)
    else:
        # Mogrify requires writing all changes to a temporary file by default,
        # but convert uses a destination file instead if one was specified. Get
        # sequences from an iterator that has generator functions wrapping it.
        # After creation, it is then copied back over the original file if all
        # tasks finish up without an exception being thrown.  This avoids
        # loading the entire sequence file up into memory.
        logging.info("Applying transformations, writing to %s",
                destination_file)
        SeqIO.write(records, destination_file, destination_file_type)


def module_function(string):
    """
    Load a function from a python module using a file name, function name
    specification of format:
        /path/to/x.py:function_name[:parameter]
    """
    parts = string.split(':', 2)
    if len(parts) < 2:
        raise ValueError(
            "Illegal specification. Should be module:function[:parameter]")
    module_path, function_name = parts[:2]

    # Import the module
    module_vars = {}
    execfile(module_path, module_vars)

    try:
        function = module_vars[function_name]
    except KeyError:
        raise argparse.ArgumentTypeError("{0} has no attribute '{1}'".format(
            module_path, function_name))

    if len(parts) == 3:
        old_function = function
        function = lambda r: old_function(r, parts[2])

    return function


def action(arguments):
    with arguments.source_file as src, \
            common.atomic_write(arguments.dest_file) as dest:
        transform_file(src, dest, arguments)
