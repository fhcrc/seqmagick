"""
Calculates the molecular weight and theoretical isoelectric
point (pI) for a set of protein sequences.
"""
import argparse
import sys

from Bio import SeqIO
from Bio.SeqUtils import ProtParam

from seqmagick import fileformat
from seqmagick import transform
from . import common


average_residue_mass = 112.5
stop_symbol = '*'
unknown_residue_char = 'X'


def remove_stop_codon(seq, stop_symbol=stop_symbol):
    if seq.endswith(stop_symbol):
        return seq[0:-1]
    return seq


def is_valid_sequence(seq, seqid):
    for i, aa in enumerate(seq):
        if aa not in ProtParam.IUPACData.protein_weights.keys():
            raise ValueError(
                "Error: Unknown residue '%s' at position %d in %s" % (
                aa, i + 1, seqid))
            #return False
    return True


def replace_unknown_residues(seq, char=unknown_residue_char):
    newseq = str(seq)
    for aa in seq:
        if aa not in ProtParam.IUPACData.protein_weights.keys():
            newseq = newseq.replace(aa, char)
    return newseq


class ProtParamCalculator(object):
    """
    Calculates the molecular weight and isoelectric point of a set of protein
    sequences.

    If allow_unknown_residues is True, residues of unknown mass
    and properties will be given an average mass of {0} Da and treated as
    untitratable for pI calculation. Otherwise, unknown residues raise a
    ValueError. sort_on can be set to 'mass' or 'pi' to sort the results on
    these parameters, in ascending order when sort_ascending=True.

    Wraps Bio.SeqUtils.ProtParam.
    """.format(average_residue_mass)

    def __init__(self, sequences, allow_unknown_residues=False, sort_on="",
                 sort_ascending=True):
        self.sequences = sequences
        self.allow_unknown_residues = allow_unknown_residues
        self.sort_on = sort_on
        self.sort_ascending = sort_ascending
        self.stats = self.calculate(sequences)

    def calculate(self, sequences):
        """
        For a list of sequences (Bio.SeqRecord objects), calculate the
        molecular weight (Mw) and isoelectric point (pI). Return a list of
        tuples (sequence, Mw, pI).
        """
        # we add a key for our unknown residue to the masses dictionary,
        # but delete it after we are done to prevent unintended side effects
        if self.allow_unknown_residues:
            ProtParam.IUPACData.protein_weights[unknown_residue_char] \
                = average_residue_mass

        stats = []
        for s in sequences:
            seq = self._prepare_sequence(s)

            params = ProtParam.ProteinAnalysis(seq)

            stats.append((s, params.molecular_weight(),
                          params.isoelectric_point()))

        # now revert this by deleting the key we added, since
        # subsequent class might rely on the absence of unknown_residue_char
        if self.allow_unknown_residues:
            ProtParam.IUPACData.protein_weights.pop(unknown_residue_char, None)

        sort_index = None
        if self.sort_on == 'mass':
            sort_index = 1
        if self.sort_on == 'pi':
            sort_index = 2
        if sort_index:
            stats = sorted(stats, key=lambda stats: stats[sort_index],
                           reverse=(not self.sort_ascending))

        self.stats = stats
        return stats

    def _prepare_sequence(self, s, mask_unknown_residues=None):
        """
        Takes a Bio.SeqRecord and returns a string version prepared for
        protparam. Removes terminal stop codons (usually '*').
        If residues of unknown mass are encountered, raises ValueError, unless
        mask_unknown_residues is True (in which case they are masked
        with unknown_residue_char, usually 'X').
        """

        if mask_unknown_residues == None:
            mask_unknown_residues = self.allow_unknown_residues

        seq = str(s.seq)
        seq = remove_stop_codon(seq)

        if mask_unknown_residues:
            seq = replace_unknown_residues(seq, char=unknown_residue_char)

        is_valid_sequence(seq, s.id)

        return seq


def build_parser(parser):
    parser.add_argument('sequence_file', help="Sequence file",
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--output-file', help="Destination trimmed file",
                        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--source-format', default=None)
    parser.add_argument('-d', '--include-description', action='store_true',
                        default=False, help="""Include the sequence description in output
            [default: %(default)s]""")
    parser.add_argument('--sort', dest='sort',
                        choices=['length-asc', 'length-desc', 'name-asc',
                                 'name-desc',
                                 'mass-asc', 'mass-desc', 'pi-asc', 'pi-desc'],
                        help='Sort output based on length, name, molecular weight or pI')
    parser.add_argument('--allow-unknown-residues', action='store_true',
                        default=False,
                        help="Allow amino acids outside the standard "
                             "20 and assign them with an an average "
                             "mass of " + str(average_residue_mass) +
                             " Da [default: %(default)s]")


def action(arguments):
    common.exit_on_sigpipe()

    # Determine file format for input and output
    source_format = (arguments.source_format or
                     fileformat.from_filename(arguments.sequence_file.name))

    with arguments.sequence_file:
        sequences = SeqIO.parse(arguments.sequence_file,
                                source_format)

        # sort based on name or length
        sorters = {'length': transform.sort_length,
                   'name': transform.sort_name, }
        directions = {'asc': 1, 'desc': 0}
        sort_ascending = True
        sort_on = ""
        if arguments.sort:
            sort_on, direction = arguments.sort.split('-')
            sort_ascending = (direction != 'desc')
            if (sort_on == 'length') or (sort_on == 'name'):
                # Sorted iterator
                key, direction = arguments.sort.split('-')
                sequences = sorters[key](source_file=arguments.sequence_file,
                                         source_file_type=source_format,
                                         direction=directions[direction])

        ppcalcer = ProtParamCalculator(sequences,
                                       arguments.allow_unknown_residues,
                                       sort_on=sort_on,
                                       sort_ascending=sort_ascending)

        stats = ppcalcer.stats

        if arguments.include_description:
            out = ((s[0].description, s[1], s[2]) for s in stats)
        else:
            out = ((s[0].id, s[1], s[2]) for s in stats)

        with arguments.output_file:
            for l in out:
                print >> arguments.output_file, "%s\t%.2f\t%.2f" % (l)
