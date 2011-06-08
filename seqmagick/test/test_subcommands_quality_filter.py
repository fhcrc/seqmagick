import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqmagick.subcommands import quality_filter


class QualityFilterTestCase(unittest.TestCase):

    def setUp(self):
        self.instance = quality_filter.QualityScoreFilter()
        self.sequence = SeqRecord(Seq('ACGT'))

    def test_nowindow_fail(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 24, 25]
        instance = quality_filter.QualityScoreFilter()
        result = instance.filter_record(self.sequence)
        self.assertEqual(None, result)

    def test_nowindow_pass(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 25, 25]
        instance = quality_filter.QualityScoreFilter()
        result = instance.filter_record(self.sequence)
        self.assertEqual(self.sequence, result)

    def test_window_pass(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 25, 25]
        instance = quality_filter.QualityScoreFilter(window_size=2)
        result = instance.filter_record(self.sequence)
        self.assertEqual(str(self.sequence), str(result))

    def test_window_truncate_noseq(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 24, 25, 25]
        instance = quality_filter.QualityScoreFilter(window_size=2)
        result = instance.filter_record(self.sequence)
        self.assertEqual(0, len(result))

    def test_window_truncate_mid(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 23, 25]
        instance = quality_filter.QualityScoreFilter(window_size=2)
        result = instance.filter_record(self.sequence)
        self.assertEqual(2, len(result))
        self.assertEqual('AC', str(result.seq))


class AmbiguousBaseFilterTestCase(unittest.TestCase):
    """
    Tests for ambiguous_base_filter
    """
    def setUp(self):
        self.records = [SeqRecord(Seq('ACGT')),
                SeqRecord(Seq('NNNN')),
                SeqRecord(Seq('NACT')),
                SeqRecord(Seq('ACGTN')),
                SeqRecord(Seq('GGNTTACT')),
                ]

    def test_drop(self):
        """
        Test that the first record (with no Ns) does not get filtered
        """
        instance = quality_filter.AmbiguousBaseFilter('drop')
        actual = list(instance.filter_records(self.records))
        self.assertEqual(1, len(actual))
        self.assertEqual(1, instance.passed)
        self.assertEqual(4, instance.failed)
        self.assertEqual(self.records[0], actual[0])

    def test_truncate(self):
        instance = quality_filter.AmbiguousBaseFilter('truncate')
        actual = list(instance.filter_records(self.records))
        self.assertEqual(5, len(actual))
        self.assertEqual(0, instance.failed)
        self.assertEqual(5, instance.passed)
        self.assertEqual(['ACGT', '', '', 'ACGT', 'GG'],
                [str(s.seq) for s in actual])

    def test_invalid_action(self):
        self.assertRaises(ValueError, quality_filter.AmbiguousBaseFilter,
                'other')
