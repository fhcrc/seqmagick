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
