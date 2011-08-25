"""
Tests for primer trim
"""
import unittest

from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqmagick.subcommands import primer_trim

class PrimerAlignerTestCase(unittest.TestCase):

    def setUp(self):
        self.primer = 'AACTGCATTTGAATGG'
        self.instance = primer_trim.PrimerAligner(self.primer, match=5.0,
                gap_open=-10.0)

    def test_max_score(self):
        self.assertEqual(len(self.primer) * 5.0, self.instance.max_score)

    def test_align_exact(self):
        sequence = ('ACTCTGTGTCACTTTAAACTGCATTTGAATGGAAGAGTAATAGTAGCAATAACGGCA'
                    'CTGATCAG')
        hamming_distance, start, end = self.instance.align(sequence)
        self.assertEqual(0, hamming_distance)
        self.assertEqual(16, start)
        self.assertEqual(31, end)

    def test_align_gap(self):
        sequence = ('ACTCTGTGTCACTTTAAACTGCATTGAATGGAAGAGTAATAGTAGCAATAACGGCA'
                    'CTGATCAG')
        hamming_distance, start, end = self.instance.align(sequence)
        expected_distance = 1
        self.assertEqual(expected_distance, hamming_distance)
        self.assertEqual(16, start)
        self.assertEqual(30, end)

class HammingDistanceTestCase(unittest.TestCase):

    def test_unequal_length(self):
        s1 = 'test'
        s2 = 'te'
        self.assertRaises(ValueError, primer_trim.hamming_distance, s1, s2)

    def test_no_difference(self):
        s1 = s2 = 'test'
        self.assertEqual(0, primer_trim.hamming_distance(s1, s2))

    def test_all_different(self):
        s1 = 'test'
        s2 = 'ACGT'
        self.assertEqual(4, primer_trim.hamming_distance(s1, s2))

    def test_basic(self):
        s1 = 'ACGT'
        s2 = 'AGGT'
        self.assertEqual(1, primer_trim.hamming_distance(s1, s2))

    def test_ambiguous(self):
        s1 = 'ACYT'
        s2 = 'ACCT'
        self.assertEqual(0, primer_trim.hamming_distance(s1, s2,
            primer_trim._iupac_ambiguous_equal))
        s2 = 'ACTT'
        self.assertEqual(0, primer_trim.hamming_distance(s1, s2,
            primer_trim._iupac_ambiguous_equal))

def _alignment_record(sequence):
    return SeqRecord(Seq(sequence,
        alphabet=Alphabet.Gapped(Alphabet.generic_dna)))

class LocatePrimersTestCase(unittest.TestCase):
    """
    Test for locate primers
    """

    def setUp(self):
        self.sequences = [_alignment_record('--A--ACTGGACGTATTC-CCCC')]

    def test_basic(self):
        forward = 'TGG'
        reverse = 'TTC'

        forward_idx, reverse_idx = primer_trim.locate_primers(self.sequences,
                forward, reverse, False, 1)

        self.assertEqual((7, 9), forward_idx)
        self.assertEqual((15, 17), reverse_idx)

    def test_no_forward(self):
        forward='GGGGGG'
        reverse = 'TTC'
        self.assertRaises(primer_trim.PrimerNotFound,
                primer_trim.locate_primers, self.sequences, forward, reverse,
                False, 1)

    def test_no_reverse(self):
        forward='TGG'
        reverse = 'GGGG'
        self.assertRaises(primer_trim.PrimerNotFound,
                primer_trim.locate_primers, self.sequences, forward, reverse,
                False, 1)

    def test_bad_order(self):
        """
        Should fail if reverse primer occurs before forward primer
        """
        reverse = 'TGG'
        forward = 'TTC'

        self.assertRaises(primer_trim.PrimerOrderError,
                primer_trim.locate_primers, self.sequences,
                forward, reverse, False, 1)
