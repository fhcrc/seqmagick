"""
Find a primer sequence in a gapped alignment, trim to amplicon
"""
import argparse
import itertools

from Bio import Alphabet, SeqIO, pairwise2
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


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
    parser.add_argument('--max-hamming-distance', type=positive(int),
            default=1, help="""Maximum Hamming distance between primer and
            alignment site (default: %(default)s)""")


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


def hamming_distance(s1, s2):
    """
    Returns the hamming distance between two strings.
    """
    if not len(s1) == len(s2):
        raise ValueError("String lengths are not equal")

    # Number of non-matching characters:
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


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

    def align(self, sequence):
        """
        Aligns the primer to the given query sequence, returning a tuple of:

            hamming_distance, start, end

        Where hamming distance is the distance between the primer and aligned
        sequence, and start and end give the start and end index of the primer
        relative to the input sequence.
        """
        seq_aln, primer_aln, score, start, end = \
                pairwise2.align.globalms(str(sequence), str(self.primer),
                        self.match, self.difference, self.gap_open,
                        self.gap_extend, one_alignment_only=True)[0]

        # Get an ungapped mapping on the sequence
        index_map = gap_index_map(seq_aln)
        ungap_map = ungap_index_map(primer_aln)

        # Trim to primer
        start = ungap_map[0]
        end = ungap_map[len(self.primer) - 1]

        ham_dist = hamming_distance(seq_aln[start:end+1],
                primer_aln[start:end+1])
        #assert primer_aln[start:end].replace('-', '') == str(self.primer)

        # TODO: handle start or end being gap
        return ham_dist, index_map[start], index_map[end]

    @property
    def max_score(self):
        """
        Maximum possible alignment score
        """
        return len(self.primer) * self.match


# Types for argparse
def iupac_ambiguous_sequence(string):
    return Seq(string, IUPAC.ambiguous_dna)


def positive(target_type):
    def inner(string):
        value = target_type(string)
        if not value >= 0:
            raise argparse.ArgumentTypeError("Invalid positive number")


def locate_primers(sequences, forward_primer, reverse_primer,
        reverse_complement, max_hamming_distance):
    """
    Find forward and reverse primers in a set of sequences, return two tuples:
    (forward_start, forward_end), (reverse_start, reverse_end)
    """
    forward_loc = None
    reverse_loc = None

    # Reverse complement the reverse primer, if appropriate
    if reverse_complement:
        reverse_primer = reverse_primer.reverse_complement()

    forward_aligner = PrimerAligner(forward_primer)
    reverse_aligner = PrimerAligner(reverse_primer)

    for sequence in sequences:
        index_map = ungap_index_map(sequence.seq)
        if forward_loc is None:
            ham_dist, start, end = \
                    forward_aligner.align(sequence.seq.ungap())
            if ham_dist <= max_hamming_distance:
                forward_loc = index_map[start], index_map[end]
        if reverse_loc is None:
            ham_dist, start, end = \
                    reverse_aligner.align(sequence.seq.ungap())
            if ham_dist <= max_hamming_distance:
                reverse_loc = index_map[start], index_map[end]
        if forward_loc and reverse_loc:
            # Both found:
            return forward_loc, reverse_loc

    # missed the forward or reverse primer:
    if not forward_loc:
        raise PrimerNotFound(forward_primer)
    else:
        raise PrimerNotFound(reverse_primer)


def action(arguments):
    """
    Trim the alignment as specified
    """
    # Load the alignment
    with arguments.source_file:
        sequences = SeqIO.parse(arguments.source_file,
                arguments.alignment_format,
                alphabet=Alphabet.Gapped(Alphabet.single_letter_alphabet))

        # Locate primers
        (forward_start, forward_end), (reverse_start, reverse_end) = \
                locate_primers(sequences, arguments.forward_primer,
                        arguments.reverse_primer, arguments.reverse_complement,
                        arguments.max_hamming_distance)

        # Generate a slice
        if arguments.include_primers:
            start = forward_start
            end = reverse_end + 1
        else:
            start = forward_end + 1
            end = reverse_start

        # Rewind
        arguments.source_file.seek(0)
        sequences = SeqIO.parse(arguments.source_file,
                arguments.alignment_format,
                alphabet=Alphabet.Gapped(Alphabet.single_letter_alphabet))

        sequences = (i[start:end] for i in sequences)
        with arguments.output_file:
            SeqIO.write(sequences, arguments.output_file,
                    arguments.alignment_format)
