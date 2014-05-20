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
            raise ValueError("Error: Unknown residue '%s' at position %d in %s" % (aa, i + 1, seqid))
            #return False
    return True


def replace_unknown_residues(seq, char=unknown_residue_char):
    newseq = str(seq)
    for aa in seq:
        if aa not in ProtParam.IUPACData.protein_weights.keys():
            newseq = newseq.replace(aa, char)
    return newseq


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
                        choices=['length-asc', 'length-desc', 'name-asc', 'name-desc',
                                 'mass-asc', 'mass-desc', 'pi-asc', 'pi-desc'],
                        help='Sort output based on length, name, molecular weight or pI')
    parser.add_argument('--allow-unknown-residues', action='store_true',
                        default=False, help="Allow amino acids outside the standard "
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
        if arguments.sort:
            sort_on, direction = arguments.sort.split('-')
            reverse = (direction == 'desc')
            if (sort_on == 'length') or (sort_on == 'name'):
                # Sorted iterator
                key, direction = arguments.sort.split('-')
                sequences = sorters[key](source_file=arguments.sequence_file,
                                         source_file_type=source_format,
                                         direction=directions[direction])

        stats = []
        for s in sequences:
            seq = str(s.seq)
            seq = remove_stop_codon(seq)

            if arguments.allow_unknown_residues:
                seq = replace_unknown_residues(seq, char=unknown_residue_char)
                ProtParam.IUPACData.protein_weights[unknown_residue_char] = average_residue_mass

            if not is_valid_sequence(seq, s.id):
                break

            params = ProtParam.ProteinAnalysis(seq)

            stats.append((s, params.molecular_weight(),
                          params.isoelectric_point()))

            if arguments.sort and sort_on == 'mass':
                stats = sorted(stats, key=lambda stats: stats[1], reverse=reverse)
            elif arguments.sort and sort_on == 'pi':
                stats = sorted(stats, key=lambda stats: stats[2], reverse=reverse)

        if arguments.include_description:
            out = ((s[0].description, s[1], s[2]) for s in stats)
        else:
            out = ((s[0].id, s[1], s[2]) for s in stats)

        with arguments.output_file:
            for l in out:
                print >> arguments.output_file, "%s\t%.2f\t%.2f" % (l)
