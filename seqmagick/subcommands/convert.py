"""
Convert between sequence formats
"""
import argparse
import logging
import os
import os.path

from Bio import SeqIO
from Bio.SeqIO import FastaIO

from seqmagick import transform
from seqmagick.fileformat import lookup_file_type

import common


def add_options(parser):
    """
    Add optional arguments to the parser
    """
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
    seq_mods.add_argument('--cut', dest='cut', metavar="start:end",
        type=common.cut_range, help='1-indexed start and end positions for '
        'cutting sequences, : separated.  Includes last item.')
    seq_mods.add_argument('--dash-gap', action='store_true', dest='dash_gap',
        help='Change . and : into - for all sequences')
    seq_mods.add_argument('--lower', action='store_true', dest='lower',
        help='Translate the sequences to lower case')
    seq_mods.add_argument('--reverse', action='store_true', dest='reverse',
        help='Reverse the order of sites in sequences')
    seq_mods.add_argument('--reverse-complement', dest='reverse_complement',
        action='store_true', help='Convert sequences into reverse complements')
    seq_mods.add_argument('--squeeze', action='store_true', dest='squeeze',
        help='Remove any gaps that are present in the same position '
        'across all sequences in an alignment')
    seq_mods.add_argument('--transcribe', dest='transcribe',
        choices=['dna2rna', 'rna2dna'],
        help='Transcription and back transcription for generic DNA and '
        'RNA. Source sequences must be the correct alphabet or this '
        'action will likely produce incorrect results.')
    seq_mods.add_argument('--translate', dest='translate',
        choices=['dna2protein', 'rna2protein',
                 'dna2proteinstop', 'rna2proteinstop'],
        help='Translate from generic DNA/RNA to proteins. Options with '
        '"stop" suffix will NOT translate through stop codons .'
        'Source sequences must be the correct alphabet or this action '
        'will likely produce incorrect results.')
    seq_mods.add_argument('--ungap', action='store_true', dest='ungap',
        help='Remove gaps in the sequence alignment')
    seq_mods.add_argument('--upper', action='store_true', dest='upper',
        help='Translate the sequences to upper case')

    seq_select = parser.add_argument_group("Record Selection")

    seq_select.add_argument('--deduplicate-sequences',
        action='store_const', const=None, default=False,
         dest='deduplicate_sequences', help='Remove any duplicate sequences '
         'by sequence content, keep the first instance seen')
    seq_select.add_argument('--deduplicated-sequences-file', action='store',
        metavar='FILE', dest='deduplicate_sequences', default=False,
        type=argparse.FileType('w'),
        help='Write all of the deduplicated sequences to a file')

    seq_select.add_argument('--deduplicate-taxa', action='store_true',
        dest='deduplicate_taxa',
        help='Remove any duplicate sequences by ID, keep the first '
        'instance seen')

    seq_select.add_argument('--head', metavar='N', dest='head', type=int,
        help='Trim down to top N sequences')
    seq_select.add_argument('--max-length', dest='max_length', metavar='N',
        type=int, help='Discard any sequences beyond the specified '
        'maximum length.  This operation occurs *before* all '
        'length-changing options such as cut and squeeze.')
    seq_select.add_argument('--min-length', dest='min_length', metavar='N',
        type=int, help='Discard any sequences less than the specified '
        'minimum length.  This operation occurs *before* all '
        'length-changing options such as cut and squeeze.')
    seq_select.add_argument('--pattern-include', metavar='regex',
        dest='pattern_include', help='Filter the sequences by '
        'regular expression in name')
    seq_select.add_argument('--pattern-exclude', metavar='regex',
        dest='pattern_exclude', help='Filter out sequences by regular '
        'expression in name')
    seq_select.add_argument('--prune-empty', action="store_true", default=False,
                        help="Prune sequences containing only gaps ('-')")
    seq_select.add_argument('--seq-pattern-include', metavar='regex',
        dest='seq_pattern_include', help='Filter the sequences by '
        'regular expression in sequence')
    seq_select.add_argument('--seq-pattern-exclude', metavar='regex',
        dest='seq_pattern_exclude', help='Filter out sequences by regular '
        'expression in sequence')
    seq_select.add_argument('--tail', metavar='N', dest='tail', type=int,
        help='Trim down to bottom N sequences')

    id_mods = parser.add_argument_group("Sequence ID Modification")
    id_mods.add_argument('--first-name', action='store_true',
        dest='first_name',
        help='Take only the first whitespace-delimited word as the '
        'name of the sequence')
    id_mods.add_argument('--name-suffix', metavar='SUFFIX',
        dest='name_suffix', help='Append a suffix to all IDs.')
    id_mods.add_argument('--name-prefix', metavar='PREFIX',
        dest='name_prefix', help='Insert a prefix for all IDs.')
    id_mods.add_argument('--pattern-replace', nargs=2,
                        metavar=('search_pattern', 'replace_pattern'),
                        help='Replace regex pattern "search_pattern" with '
                             '"replace_pattern" in sequence ID')
    id_mods.add_argument('--strip-range', dest='strip_range',
        action='store_true', help='Strip ranges from sequences IDs, '
        'matching </x-y>')

    format_group = parser.add_argument_group('Format Options')
    format_group.add_argument('--input-format', metavar='Format',
            help="Input file format (default: determine from extension)")
    format_group.add_argument('--output-format', metavar='Format',
            help="Output file format (default: determine from extension)")

    return parser


def build_parser(parser):
    """
    Add shared arguments to the convert or mogrify parser.
    """
    add_options(parser)
    parser.add_argument('source_file', type=common.sequence_file,
            help="Input sequence file")
    parser.add_argument('dest_file', help="Output file")

    return parser


def transform_file(source_file, destination_file, arguments):

    # Get just the file name, useful for naming the temporary file.
    file_ext = os.path.splitext(source_file)[1]
    source_file_type = (arguments.input_format or lookup_file_type(file_ext))

    output_ext = os.path.splitext(destination_file)[1]
    destination_file_type = (arguments.output_format or
            lookup_file_type(output_ext))

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
        records = SeqIO.parse(source_file, source_file_type)


    #########################################
    # Apply generator functions to iterator.#
    #########################################

    logging.info('Setting up transform functions for file: %s', source_file)

    # Deduplication occurs first, to get a checksum of the
    # original sequence and to store the id field before any
    # transformations occur.

    if arguments.max_length:
        records = transform.max_length_discard(records, arguments.max_length)

    if arguments.min_length:
        records = transform.min_length_discard(records, arguments.min_length)

    if (arguments.deduplicate_sequences or
            arguments.deduplicate_sequences is None):
        records = transform.deduplicate_sequences(
            records, arguments.deduplicate_sequences)

    if arguments.deduplicate_taxa:
        records = transform.deduplicate_taxa(records)

    if arguments.dash_gap:
        records = transform.dashes_cleanup(records)

    if arguments.first_name:
        records = transform.first_name_capture(records)
    if arguments.upper:
        records = transform.upper_sequences(records)

    if arguments.lower:
        records = transform.lower_sequences(records)

    if arguments.prune_empty:
        records = transform.prune_empty(records)

    if arguments.reverse:
        records = transform.reverse_sequences(records)

    if arguments.reverse_complement:
        records = transform.reverse_complement_sequences(records)

    if arguments.ungap:
        records = transform.ungap_sequences(records)

    if arguments.name_prefix:
        records = transform.name_insert_prefix(records, arguments.name_prefix)

    if arguments.name_suffix:
        records = transform.name_append_suffix(records, arguments.name_suffix)

    if arguments.pattern_include:
        records = transform.name_include(records, arguments.pattern_include)

    if arguments.pattern_exclude:
        records = transform.name_exclude(records, arguments.pattern_exclude)

    if arguments.pattern_replace:
        search_pattern, replace_pattern = arguments.pattern_replace
        records = transform.name_replace(records, search_pattern,
                replace_pattern)

    if arguments.head and arguments.tail:
        raise ValueError("Error: head and tail are mutually exclusive "
                "at the moment.")

    if arguments.head:
        records = transform.head(records, arguments.head)

    if arguments.strip_range:
        records = transform.strip_range(records)

    if arguments.tail:
        # To know where to begin including records for tail, we need to count
        # the total number of records, which requires going through the entire
        # file and additional time.
        record_count = sum(1 for record in SeqIO.parse(source_file, source_file_type))
        records = transform.tail(records, arguments.tail, record_count)

    if arguments.transcribe:
        records = transform.transcribe(records, arguments.transcribe)

    if arguments.translate:
        records = transform.translate(records, arguments.translate)

    if arguments.squeeze:
        logging.info("Performing squeeze")
        gaps = []
        # Need to iterate an additional time to determine which
        # gaps are share between all sequences in an alignment.
        for record in SeqIO.parse(source_file, source_file_type):
            # Use numpy to prepopulate a gaps list.
            if len(gaps) == 0:
                gaps_length = len(record.seq)
                #gaps = list(ones( (gaps_length), dtype=int16 ))
                gaps = [1] * gaps_length
            gaps = map(transform.gap_check, gaps, list(str(record.seq)))
        records = transform.squeeze(records, gaps)
        logging.debug("Squeeze gaps:\n%s", gaps)

    # cut needs to go after squeeze or the gaps list will no longer be relevent.
    # It is probably best not to use squeeze and cut together in most cases.
    if arguments.cut:
        start, end = arguments.cut
        records = transform.cut_sequences(records, start=start, end=end)

   # Only the fasta format is supported, as SeqIO.write does not have a 'wrap' parameter.
    if (arguments.line_wrap is not None and destination_file_type == 'fasta'
            and source_file_type == 'fasta'):
        logging.info("Attempting to write fasta with %d line breaks.",
                arguments.line_wrap)

        with open(destination_file, "w") as handle:
            writer = FastaIO.FastaWriter(handle, wrap=arguments.line_wrap)
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


def action(arguments):
    transform_file(arguments.source_file, arguments.dest_file, arguments)
