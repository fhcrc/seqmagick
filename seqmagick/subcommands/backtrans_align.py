"""
Given a protein alignment and unaligned nucleotides, align the nucleotides
using the protein alignment.

Protein and nucleotide sequence files must contain the same number of
sequences, in the same order, with the same IDs.
"""

# TODO: Add tests

import itertools
import logging
import sys

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqmagick import fileformat

from . import common

TRANSLATION_TABLES = {
    'standard': CodonTable.unambiguous_dna_by_name["Standard"],
    'standard-ambiguous': CodonTable.ambiguous_dna_by_name["Standard"],
    'vertebrate-mito': CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
}


def build_parser(parser):
    parser.add_argument(
        'protein_align', type=common.FileType('r'), help='Protein Alignment')
    parser.add_argument(
        'nucl_align', type=common.FileType('r'), help='FASTA Alignment')
    parser.add_argument(
        '-o', '--out-file', type=common.FileType('w'),
        default=sys.stdout, metavar='destination_file',
        help='Output destination. Default: STDOUT')
    parser.add_argument(
        '-t', '--translation-table', choices=TRANSLATION_TABLES,
        default='standard-ambiguous',
        help='Translation table to use. [Default: %(default)s]')
    parser.add_argument(
        '-a', '--fail-action', choices=('fail', 'warn', 'none'), default='fail',
        help='Action to take on an ambiguous codon [default: %(default)s]')

    return parser


def batch(iterable, chunk_size):
    """
    Return items from iterable in chunk_size bits.

    If len(iterable) % chunk_size > 0, the last item returned will be shorter.
    """
    i = iter(iterable)
    while True:
        r = list(itertools.islice(i, chunk_size))
        if not r:
            break
        yield r


class AlignmentMapper(object):
    def __init__(self, translation_table, unknown_action='fail'):
        self.translation_table = translation_table
        self.unknown_action = unknown_action

    def _validate_translation(self, aligned_prot, aligned_nucl):
        """
        Given a seq for protein and nucleotide, ensure that the translation holds
        """
        codons = [''.join(i) for i in batch(str(aligned_nucl), 3)]
        for codon, aa in zip(codons, str(aligned_prot)):
            # Check gaps
            if codon == '---' and aa == '-':
                continue

            try:
                trans = self.translation_table.forward_table[codon]
                if not trans == aa:
                    raise ValueError("Codon {0} translates to {1}, not {2}".format(
                        codon, trans, aa))
            except (KeyError, CodonTable.TranslationError):
                if aa != 'X':
                    if self.unknown_action == 'fail':
                        raise ValueError("Unknown codon: {0} mapped to {1}".format(
                            codon, aa))
                    elif self.unknown_action == 'warn':
                        logging.warn('Cannot verify that unknown codon %s '
                                     'maps to %s', codon, aa)
        return True

    def map_alignment(self, prot_seq, nucl_seq):
        """
        Use aligned prot_seq to align nucl_seq
        """
        if prot_seq.id != nucl_seq.id:
            logging.warning(
                'ID mismatch: %s != %s. Are the sequences in the same order?',
                prot_seq.id, nucl_seq.id)

        # Ungap nucleotides
        codons = batch(str(nucl_seq.seq.replace('-', '')), 3)
        codons = [''.join(i) for i in codons]
        codon_iter = iter(codons)

        ungapped_prot = str(prot_seq.seq).replace('-', '')

        if len(ungapped_prot) != len(codons):
            table = self.translation_table.forward_table
            prot_str = ' '.join(' ' + p + ' ' for p in ungapped_prot)
            codon_str = ' '.join(codons)
            trans_str = ' '.join(' ' + table.get(codon, 'X') + ' '
                                 for codon in codons)
            raise ValueError("""Length of codon sequence ({0}) does not match \
length of protein sequence ({1}) for {2}
Protein:       {3}
Codons:        {4}
Trans. Codons: {5}""".format(len(codons), len(ungapped_prot), nucl_seq.id, prot_str,
                            codon_str, trans_str))

        try:
            nucl_align = ['---' if p == '-' else next(codon_iter)
                          for p in str(prot_seq.seq)]
        except StopIteration:
            assert False  # Should be checked above

        result = SeqRecord(Seq(''.join(nucl_align)), id=nucl_seq.id,
                           description=nucl_seq.description)

        # Validate
        self._validate_translation(prot_seq.seq.upper(), result.seq.upper())

        return result

    def map_all(self, prot_alignment, nucl_sequences):
        """
        Convert protein sequences to nucleotide alignment
        """
        zipped = itertools.zip_longest(prot_alignment, nucl_sequences)
        for p, n in zipped:
            if p is None:
                raise ValueError("Exhausted protein sequences")
            elif n is None:
                raise ValueError("Exhausted nucleotide sequences")
            yield self.map_alignment(p, n)

def action(arguments):
    """
    Run
    """
    # Ignore SIGPIPE, for head support
    common.exit_on_sigpipe()
    logging.basicConfig()

    prot_sequences = SeqIO.parse(arguments.protein_align,
                                 fileformat.from_handle(arguments.protein_align))
    nucl_sequences = SeqIO.parse(arguments.nucl_align,
                                 fileformat.from_handle(arguments.nucl_align))

    instance = AlignmentMapper(TRANSLATION_TABLES[arguments.translation_table],
                               arguments.fail_action)

    SeqIO.write(instance.map_all(prot_sequences, nucl_sequences),
                arguments.out_file, fileformat.from_filename(arguments.out_file.name))
