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
        self.assertEquals(None, result)

    def test_nowindow_pass(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 25, 25]
        instance = quality_filter.QualityScoreFilter()
        result = instance.filter_record(self.sequence)
        self.assertEquals(self.sequence, result)

    def test_window_pass(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 25, 25]
        instance = quality_filter.QualityScoreFilter(window_size=2)
        result = instance.filter_record(self.sequence)
        self.assertEquals(str(self.sequence), str(result))

    def test_window_truncate_noseq(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 24, 25, 25]
        instance = quality_filter.QualityScoreFilter(window_size=2)
        result = instance.filter_record(self.sequence)
        self.assertEquals(0, len(result))

    def test_window_truncate_mid(self):
        self.sequence.letter_annotations['phred_quality'] = [25, 25, 23, 25]
        instance = quality_filter.QualityScoreFilter(window_size=2)
        result = instance.filter_record(self.sequence)
        self.assertEquals(2, len(result))
        self.assertEquals('AC', str(result.seq))

