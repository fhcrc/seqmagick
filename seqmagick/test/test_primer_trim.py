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
        self.assertEquals(len(self.primer) * 5.0, self.instance.max_score)

    def test_align_exact(self):
        sequence = ('ACTCTGTGTCACTTTAAACTGCATTTGAATGGAAGAGTAATAGTAGCAATAACGGCA'
                    'CTGATCAG')
        hamming_distance, start, end = self.instance.align(sequence)
        self.assertEquals(0, hamming_distance)
        self.assertEquals(16, start)
        self.assertEquals(31, end)

    def test_align_gap(self):
        sequence = ('ACTCTGTGTCACTTTAAACTGCATTGAATGGAAGAGTAATAGTAGCAATAACGGCA'
                    'CTGATCAG')
        hamming_distance, start, end = self.instance.align(sequence)
        expected_distance = 1
        self.assertEquals(expected_distance, hamming_distance)
        self.assertEquals(16, start)
        self.assertEquals(30, end)


class HammingDistanceTestCase(unittest.TestCase):

    def test_unequal_length(self):
        s1 = 'test'
        s2 = 'te'
        self.assertRaises(ValueError, primer_trim.hamming_distance, s1, s2)

    def test_no_difference(self):
        s1 = s2 = 'test'
        self.assertEquals(0, primer_trim.hamming_distance(s1, s2))

    def test_all_different(self):
        s1 = 'test'
        s2 = 'ACGT'
        self.assertEquals(4, primer_trim.hamming_distance(s1, s2))

    def test_basic(self):
        s1 = 'ACGT'
        s2 = 'AGGT'
        self.assertEquals(1, primer_trim.hamming_distance(s1, s2))

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

        self.assertEquals((7, 9), forward_idx)
        self.assertEquals((15, 17), reverse_idx)

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


class IsolateRegionTestCase(unittest.TestCase):

    def setUp(self):
        self.sequences = [_alignment_record('--A--ACTGGACGTATTC-CCCC'),
                          _alignment_record('--AGCACTGGA---ATTC-CCCC')]

    def test_no_isolation(self):
        result = list(primer_trim.isolate_region(self.sequences, 0,
            len(self.sequences[0])))

        self.assertEquals(self.sequences, result)

    def test_single_loc(self):
        start = 2
        end = 3
        result = list(primer_trim.isolate_region(self.sequences, start, end))
        for seq in result:
            self.assertEqual('--A--------------------', str(seq.seq))

    def test_middle(self):
        expected = ['--A--ACTGGA------------', '--AGCACTGGA------------']
        start = 1
        end = 11

        actual = list(primer_trim.isolate_region(self.sequences, start, end))
        actual = [str(s.seq) for s in actual]
        self.assertEquals(expected, actual)

    def test_invalid(self):
        self.assertRaises(ValueError, primer_trim.isolate_region(
                self.sequences, 5, 5).next)
        self.assertRaises(ValueError, primer_trim.isolate_region(
                self.sequences, 10, 5).next)
