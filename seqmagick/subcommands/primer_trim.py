"""
Find a primer sequence in a gapped alignment, trim to amplicon
"""
import argparse
import itertools
import logging
import operator
import sys

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq

from seqmagick import transform, fileformat

from . import common


def build_parser(parser):
    parser.add_argument(
        'source_file', type=argparse.FileType('r'), help="Source alignment file")
    parser.add_argument(
        'output_file', type=argparse.FileType('w'), help="Destination trimmed file")
    parser.add_argument(
        'forward_primer', type=iupac_ambiguous_sequence,
        help="The forward primer used")
    parser.add_argument(
        'reverse_primer', type=iupac_ambiguous_sequence,
        help="""The reverse primer used. By default the reverse primer
            is assumed to be a subsequence of the top strand (that is,
            the reverse complement of an actual downstream PCR
            primer). Use --reverse-is-revcomp if this is not the
            case.""")
    parser.add_argument(
        '--reverse-is-revcomp', dest="reverse_complement", default=False,
        action='store_true', help="""Reverse primer is written as the
            reverse complement of the top strand (default:
            %(default)s)""")
    parser.add_argument(
        '--source-format', default=None,
        help='Alignment format (default: detect from extension')
    parser.add_argument(
        '--output-format', default=None,
        help='Alignment format (default: detect from extension')
    parser.add_argument(
        '--include-primers', action="store_true", default=False,
        help="""Include the primers in the output (default: %(default)s)""")
    parser.add_argument(
        '--max-hamming-distance', type=common.positive_value(int), default=1,
        help="""Maximum Hamming
            distance between primer and alignment site (default: %(default)s).
            IUPAC ambiguous bases in the primer matching unambiguous bases in
            the alignment are not penalized""")
    parser.add_argument(
        '--prune-action', default='trim', choices=list(_ACTIONS.keys()),
        help="""Action to take. Options are trim (trim to the region
            defined by the two primers, decreasing the width of the alignment),
            or isolate (convert all characters outside the primer-defined area
            to gaps). default: %(default)s""")


# Sequence-related functions
def ungap_index_map(sequence, gap_chars='-'):
    """
    Returns a dict mapping from an index in the ungapped sequence to an index
    in the gapped sequence.

    >>> ungap_index_map('AC-TG-')
    {0: 0, 1: 1, 2: 3, 3: 4}
    """
    counter = itertools.count(0).__next__
    ungap_indexes = [
        counter() if c not in gap_chars else None for c in iter(sequence)
    ]
    return dict(
        (ungapped, gapped)
        for ungapped, gapped in zip(ungap_indexes, range(len(sequence)))
        if ungapped is not None)


def gap_index_map(sequence, gap_chars='-'):
    """
    Opposite of ungap_index_map: returns mapping from gapped index to ungapped
    index.

    >>> gap_index_map('AC-TG-')
    {0: 0, 1: 1, 3: 2, 4: 3}
    """
    return dict(
        (v, k) for k, v in list(ungap_index_map(sequence, gap_chars).items()))


def _iupac_ambiguous_equal(ambig_base, unambig_base):
    """
    Tests two bases for equality, accounting for IUPAC ambiguous DNA

    ambiguous base may be IUPAC ambiguous, unambiguous must be one of ACGT
    """
    iupac_translation = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'U': 'U',
        'R': 'AG',
        'Y': 'CT',
        'S': 'GC',
        'W': 'AT',
        'K': 'GT',
        'M': 'AC',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT',
        '-': '-'
    }
    for i in (ambig_base, unambig_base):
        if not len(i) == 1:
            raise ValueError("only one base may be passed.")

    return unambig_base.upper() in iupac_translation[ambig_base.upper()]


def hamming_distance(s1, s2, equality_function=operator.eq):
    """
    Returns the hamming distance between two strings.
    """
    if not len(s1) == len(s2):
        raise ValueError("String lengths are not equal")

    # Number of non-matching characters:
    return sum(not equality_function(c1, c2) for c1, c2 in zip(s1, s2))


class PrimerNotFound(Exception):
    pass


class PrimerOrderError(Exception):
    def __init__(self, forward_indexes, reverse_indexes):
        super(PrimerOrderError, self).__init__(
            "Reverse primer before forward primer: {0} > {1}".format(
                forward_indexes, reverse_indexes))


class PrimerAligner(object):
    """
    Get positions of pairwise alignments of a primer to a sequence.
    """

    def __init__(self, primer, match=5, difference=-4, gap_open=-10,
                 gap_extend=-0.5, penalize_end_gaps=False):
        self.primer = primer
        self.match = match
        self.difference = difference
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.penalize_end_gaps = penalize_end_gaps

    def align(self, sequence):
        """
        Aligns the primer to the given query sequence, returning a tuple of:

            hamming_distance, start, end

        Where hamming distance is the distance between the primer and aligned
        sequence, and start and end give the start and end index of the primer
        relative to the input sequence.
        """
        seq_aln, primer_aln, score, start, end = pairwise2.align.globalms(
            str(sequence).upper(), str(self.primer).upper(),
            self.match, self.difference, self.gap_open,
            self.gap_extend, one_alignment_only=True,
            penalize_end_gaps=self.penalize_end_gaps)[0]

        # Get an ungapped mapping on the sequence
        index_map = gap_index_map(seq_aln)
        ungap_map = ungap_index_map(primer_aln)

        # Trim to primer
        start = ungap_map[0]
        end = ungap_map[len(self.primer) - 1]

        trimmed = seq_aln[start:end + 1]

        ham_dist = hamming_distance(primer_aln[start:end + 1], trimmed,
                                    _iupac_ambiguous_equal)
        # assert primer_aln[start:end].replace('-', '') == str(self.primer)

        # TODO: handle start or end being gap better. For now, just give up
        # and return maxint for the hamming distance
        if trimmed.endswith('-'):
            tail = len(trimmed) - len(trimmed.rstrip('-'))
            end = index_map[end - tail] + 1
            ham_dist = sys.maxsize
        else:
            end = index_map[end]
        if trimmed.startswith('-'):
            start = 0
            ham_dist = sys.maxsize
        else:
            start = index_map[start]

        return ham_dist, start, end

    @property
    def max_score(self):
        """
        Maximum possible alignment score
        """
        return len(self.primer) * self.match


# Types for argparse
def iupac_ambiguous_sequence(string):
    return Seq(string, IUPAC.ambiguous_dna)


def locate_primers(sequences, forward_primer, reverse_primer,
                   reverse_complement, max_hamming_distance):
    """
    Find forward and reverse primers in a set of sequences, return two tuples:
    (forward_start, forward_end), (reverse_start, reverse_end)
    """
    forward_loc = None
    reverse_loc = None
    seq_length = None

    # Reverse complement the reverse primer, if appropriate
    if reverse_complement:
        reverse_primer = reverse_primer.reverse_complement()

    forward_aligner = PrimerAligner(forward_primer)
    reverse_aligner = PrimerAligner(reverse_primer)

    for i, sequence in enumerate(sequences):
        if seq_length is None:
            seq_length = len(sequence)
        elif len(sequence) != seq_length:
            raise ValueError(("Sequence Length Heterogeneity: {0} != {1}. "
                              "Is this an alignment?").format(
                                  len(sequence), seq_length))
        index_map = ungap_index_map(sequence.seq)
        if forward_loc is None:
            ham_dist, start, end = forward_aligner.align(sequence.seq.ungap())
            if ham_dist <= max_hamming_distance:
                forward_loc = index_map[start], index_map[end]
                logging.info("Forward in sequence %d: indexes %d to %d", i + 1,
                             *forward_loc)
        if reverse_loc is None:
            ham_dist, start, end = reverse_aligner.align(sequence.seq.ungap())
            if ham_dist <= max_hamming_distance:
                reverse_loc = index_map[start], index_map[end]
                logging.info("Reverse in sequence %d: indexes %d to %d", i + 1,
                             *reverse_loc)
        if forward_loc and reverse_loc:
            # Both found
            # Check order
            if forward_loc[0] > reverse_loc[0]:
                raise PrimerOrderError(forward_loc[0], reverse_loc[0])
            return forward_loc, reverse_loc
        else:
            logging.debug(
                "Sequence %d: %d/2 primers found", i + 1,
                sum(j is not None for j in (forward_loc, reverse_loc)))

    # Did not find either the forward or reverse primer:
    if not forward_loc:
        raise PrimerNotFound(forward_primer)
    else:
        raise PrimerNotFound(reverse_primer)


def trim(sequences, start, end):
    """
    Slice the input sequences from start to end
    """
    logging.info("Trimming from %d to %d", start, end)
    return (sequence[start:end] for sequence in sequences)


# Prune actions
_ACTIONS = {'trim': trim, 'isolate': transform.isolate_region}


def action(arguments):
    """
    Trim the alignment as specified
    """
    # Determine file format for input and output
    source_format = (arguments.source_format or
                     fileformat.from_handle(arguments.source_file))
    output_format = (arguments.output_format or
                     fileformat.from_handle(arguments.output_file))

    # Load the alignment
    with arguments.source_file:
        sequences = SeqIO.parse(
            arguments.source_file,
            source_format)

        # Locate primers
        (forward_start, forward_end), (reverse_start, reverse_end) = locate_primers(
            sequences, arguments.forward_primer,
            arguments.reverse_primer, arguments.reverse_complement,
            arguments.max_hamming_distance)

        # Generate slice indexes
        if arguments.include_primers:
            start = forward_start
            end = reverse_end + 1
        else:
            start = forward_end + 1
            end = reverse_start

        # Rewind the input file
        arguments.source_file.seek(0)
        sequences = SeqIO.parse(
            arguments.source_file,
            source_format)

        # Apply the transformation
        prune_action = _ACTIONS[arguments.prune_action]
        transformed_sequences = prune_action(sequences, start, end)

        with arguments.output_file:
            SeqIO.write(transformed_sequences, arguments.output_file,
                        output_format)
