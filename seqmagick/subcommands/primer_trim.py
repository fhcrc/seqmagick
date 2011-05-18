"""
Find a primer sequence in a gapped alignment, trim to amplicon
"""
import argparse
import itertools

from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2


def ungap_index_map(sequence, gap_chars='-'):
    """
    Returns a dict mapping from an index in the ungapped sequence to an index
    in the gapped sequence.

    >>> ungap_index_map('AC-TG-')
    {0: 0, 1: 1, 2: 3, 3: 4}
    """
    counter = itertools.count(0).next
    ungap_indexes = [counter() if c not in gap_chars else None
                     for c in iter(sequence)]
    return  dict((ungapped, gapped)
                 for ungapped, gapped in zip(ungap_indexes,
                                             xrange(len(sequence)))
                 if ungapped is not None)

def gap_index_map(sequence, gap_chars='-'):
    """
    Opposite of ungap_index_map: returns mapping from gapped index to ungapped
    index.

    >>> gap_index_map('AC-TG-')
    {0: 0, 1: 1, 3: 2, 4: 3}
    """
    return dict((v, k)
                for k, v in ungap_index_map(sequence, gap_chars).items())


class PrimerNotFound(Exception):
    pass


class PrimerAligner(object):
    """
    Get positions of pairwise alignments of a primer to a sequence.
    """
    def __init__(self, primer, match=5, difference=-4, gap_open=-10,
            gap_extend=-0.5):
        self.primer = primer
        self.match = match
        self.difference = difference
        self.gap_open = gap_open
        self.gap_extend = gap_extend

    def align(self, sequence, relative_score=False):
        seq_aln, primer_aln, score, start, end = \
                pairwise2.align.localms(str(sequence), str(self.primer),
                        self.match, self.difference, self.gap_open,
                        self.gap_extend, one_alignment_only=True)[0]

        # Get an ungapped mapping on the sequence
        index_map = gap_index_map(seq_aln)

        relative_score = score / float(self.max_score)

        # TODO: handle start or end being gap
        return score, relative_score, index_map[start], index_map[end]

    @property
    def max_score(self):
        """
        Maximum possible alignment score
        """
        return len(self.primer) * self.match


def iupac_ambiguous_sequence(string):
    return Seq(string, IUPAC.ambiguous_dna)

def proportion(string):
    f = float(string)
    if not 0.0 < f <= 1.0:
        raise argparse.ArgumentTypeError("Invalid proportion")


def build_parser(parser):
    parser.add_argument('source_file', help="Source alignment file",
            type=argparse.FileType('r'))
    parser.add_argument('output_file', help="Destination trimmed file",
            type=argparse.FileType('w'))
    parser.add_argument('forward_primer',
            help="The forward primer used", type=iupac_ambiguous_sequence)
    parser.add_argument('reverse_primer',
            help="The reverse primer used", type=iupac_ambiguous_sequence)
    parser.add_argument('--reverse-complement', default=False,
            action='store_true', help="""Reverse primer is 5' to 3' (requires
            reverse reverse complement.""", dest="reverse_complement")
    parser.add_argument('--alignment-format', default='fasta',
            help='Alignment format (default: %(default)s)')
    parser.add_argument('--include-primers', default=False,
            help='Include the primers in the output (default: %(default)s)')
    parser.add_argument('--min-relative-score', type=proportion, default=0.80,
            help="Minimum relative score to consider a match (default: %(default)s)")


def locate_primers(sequences, forward_primer, reverse_primer,
        reverse_complement, score_threshold):
    """
    Find forward and reverse primers in a set of sequences, return two tuples:
    (forward_start, forward_end), (reverse_start, reverse_end)
    """
    forward_loc = reverse_loc = None

    forward_aligner = PrimerAligner(forward_primer)
    reverse_aligner = PrimerAligner(reverse_primer)

    # Reverse complement the reverse primer, if appropriate
    if reverse_complement:
        reverse_primer.reverse_complement()

    for sequence in sequences:
        pass
    raise NotImplementedError()


def action(arguments):
    """

    """
    # Load the alignment
    with arguments.source_file:
        sequences = SeqIO.parse(arguments.source_file,
                arguments.alignment_format,
                alphabet=Alphabet.Gapped(Alphabet.single_letter_alphabet))

        # Locate primers
        forward_start, forward_end, reverse_start, reverse_end = \
                locate_primers(sequences, arguments.forward_primer,
                        arguments.reverse_primer, arguments.reverse_complement,
                        arguments.min_relative_score)

    raise NotImplementedError()
