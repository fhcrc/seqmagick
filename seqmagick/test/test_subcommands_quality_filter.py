import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqmagick.subcommands import quality_filter

class QualityFilterTestCase(unittest.TestCase):

    def setUp(self):
        self.instance = quality_filter.QualityScoreFilter(25.0)
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

class WindowQualityFilterTestCase(unittest.TestCase):

    def setUp(self):
        self.instance = quality_filter.WindowQualityScoreFilter(2, 25)
        self.sequence = SeqRecord(Seq('ACGT'))

    def test_window_pass(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 25, 25]
        result = self.instance.filter_record(self.sequence)
        self.assertEqual(str(self.sequence), str(result))

    def test_window_truncate_noseq(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 24, 25, 25]
        result = self.instance.filter_record(self.sequence)
        self.assertEqual(0, len(result))

    def test_window_truncate_mid(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 23, 25]
        result = self.instance.filter_record(self.sequence)
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

class MaxAmbiguousFilterTestCase(unittest.TestCase):
    """
    """
    def setUp(self):
        self.records = [SeqRecord(Seq('ACGT')),
                SeqRecord(Seq('NNNN')),
                SeqRecord(Seq('NACT')),
                SeqRecord(Seq('ACNTN')),
                SeqRecord(Seq('GGNTTNACT')),
                ]

    def test_none(self):
        instance = quality_filter.MaxAmbiguousFilter(0)
        filtered = list(instance.filter_records(self.records))
        self.assertEqual(len(filtered), 1)
        self.assertEqual(str(self.records[0].seq), str(filtered[0].seq))

    def test_10(self):
        instance = quality_filter.MaxAmbiguousFilter(10)
        filtered = list(instance.filter_records(self.records))
        self.assertEqual(filtered, self.records)

    def test_1(self):
        instance = quality_filter.MaxAmbiguousFilter(1)
        filtered = list(instance.filter_records(self.records))
        self.assertEqual([self.records[i] for i in (0, 2)], filtered)


class MinLengthFilterTestCase(unittest.TestCase):

    def setUp(self):
        self.sequences = [SeqRecord(Seq('ACGT')),
                          SeqRecord(Seq('ACTTT')), ]

    def test_none_pass(self):
        instance = quality_filter.MinLengthFilter(6)
        actual = list(instance.filter_records(self.sequences))
        self.assertEqual([], actual)

    def test_all_pass(self):
        instance = quality_filter.MinLengthFilter(4)
        actual = list(instance.filter_records(self.sequences))
        self.assertEqual(self.sequences, actual)

    def test_some_pass(self):
        instance = quality_filter.MinLengthFilter(5)
        actual = list(instance.filter_records(self.sequences))
        self.assertEqual(self.sequences[1:], actual)

class MaxLengthFilterTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = [SeqRecord(Seq('ACGT')),
                          SeqRecord(Seq('ACTTT')), ]

    def test_none_truncated(self):
        instance = quality_filter.MaxLengthFilter(6)
        actual = list(instance.filter_records(self.sequences))
        self.assertEqual(self.sequences, actual)

    def test_some_truncated(self):
        instance = quality_filter.MaxLengthFilter(4)
        actual = list(instance.filter_records(self.sequences))
        self.assertEqual(['ACGT', 'ACTT'], [str(s.seq) for s in actual])

    def test_all_truncated(self):
        instance = quality_filter.MaxLengthFilter(3)
        actual = list(instance.filter_records(self.sequences))
        self.assertEqual(['ACG', 'ACT'], [str(s.seq) for s in actual])
        self.assertEqual([i.id for i in self.sequences], [i.id for i in actual])
