"""
Tests for primer trim
"""

import unittest

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
        score, relative_score, start, end = self.instance.align(sequence)
        self.assertEquals(1.0, relative_score)
        self.assertEquals(self.instance.max_score, score)
        self.assertEquals(16, start)
        self.assertEquals(32, end)

    def test_align_gap(self):
        sequence = ('ACTCTGTGTCACTTTAAACTGCATTGAATGGAAGAGTAATAGTAGCAATAACGGCA'
                    'CTGATCAG')
        score, relative_score, start, end = self.instance.align(sequence)
        max_score = self.instance.max_score
        # Expected: maximum - 1 match - gapopen
        expected_score = max_score - 10 - 5
        self.assertEquals(expected_score / max_score, relative_score)
        self.assertEquals(expected_score, score)
        self.assertEquals(16, start)
        self.assertEquals(31, end)

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

