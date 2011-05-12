"""
Tests for seqmagick.magickwrap
"""

import StringIO
import unittest
import tempfile

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from seqmagick import magickwrap


def seqrecord(sequence_id, sequence_text, alphabet=None):
    """
    Quick shortcut to make a SeqRecord
    """
    return SeqRecord(Seq(sequence_text, alphabet), id=sequence_id)

class MagickWrapMixin(object):
    """
    Mix-in which automates creating a MagickWrap instance.
    """

    def setUp(self):
        self.infile = StringIO.StringIO()
        self.outfile = StringIO.StringIO()
        self.instance = magickwrap.MagickWrap(tempfile.gettempdir(),
                [self.infile], self.outfile)

    def tearDown(self):
        pass

class PatternReplaceTestCase(MagickWrapMixin, unittest.TestCase):

    def create_sequences(self):
        return [
    seqrecord('test_sequence_1', 'ACTGT'),
    seqrecord('test_REPLACE_2', 'ACTGT'),
    seqrecord('other_sequence', 'ATGAG'),
    ]

    def setUp(self):
        super(PatternReplaceTestCase, self).setUp()
        self.sequences = self.create_sequences()

    def tearDown(self):
        super(PatternReplaceTestCase, self).tearDown()

    def test_pattern_replace_none(self):
        result = self.instance._name_replace(self.sequences, 'ZZZ', 'MATCH')
        result = list(result)
        self.assertEqual(self.sequences, result)

    def test_pattern_replace_static(self):
        result = self.instance._name_replace(self.sequences, '_REPLACE_',
                '_DONE_')
        result = list(result)
        expected = self.create_sequences()
        expected[1].id = 'test_DONE_2'
        self.assertEqual(self.sequences, result)

    def test_pattern_replace_case_insensitive(self):
        """
        Substitutions are case insensitive
        """
        result = self.instance._name_replace(self.sequences, '_replace_',
                '_DONE_')
        result = list(result)
        expected = self.create_sequences()
        expected[1].id = 'test_DONE_2'
        self.assertEqual(self.sequences, result)

    def test_pattern_replace_group(self):
        """
        Make sure capturing groups work
        """
        result = self.instance._name_replace(self.sequences, '_(repl)ace_',
                '_DONE-\\1_')
        result = list(result)
        expected = self.create_sequences()
        expected[1].id = 'test_DONE-repl_2'
        self.assertEqual(self.sequences, result)

class SqueezeTestCase(MagickWrapMixin, unittest.TestCase):

    def setUp(self):
        super(SqueezeTestCase, self).setUp()

        self.sequences = [
            seqrecord('sequence_1', 'AC-G--'),
            seqrecord('sequence_2', '-C-GT-'),
            seqrecord('sequence_3', '-T-AG-'),
        ]

    def test_basic_squeeze(self):
        result = list(self.instance._squeeze(self.sequences,
            [False, False, True, False, False]))

        self.assertEqual([4, 4, 4], [len(i) for i in result])
        self.assertEqual([i.id for i in self.sequences], [i.id for i in result])
        expected = [
            seqrecord('sequence_1', 'ACG-'),
            seqrecord('sequence_2', '-CGT'),
            seqrecord('sequence_3', '-TAG'),
        ]

        self.assertEqual([str(i.seq) for i in expected],
                [str(i.seq) for i in result])
