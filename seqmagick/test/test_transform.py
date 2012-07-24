"""
Tests for seqmagick.transform
"""

from cStringIO import StringIO
import functools
import logging
import unittest

from Bio import Alphabet, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from seqmagick import transform

logging.basicConfig(level=logging.FATAL)

def _alignment_record(sequence):
    return SeqRecord(Seq(sequence,
        alphabet=Alphabet.Gapped(Alphabet.generic_dna)))

def seqrecord(sequence_id, sequence_text, alphabet=Alphabet.generic_dna):
    """
    Quick shortcut to make a SeqRecord
    """
    return SeqRecord(Seq(sequence_text, alphabet), id=sequence_id)

class PatternReplaceTestCase(unittest.TestCase):

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
        result = transform.name_replace(self.sequences, 'ZZZ', 'MATCH')
        result = list(result)
        self.assertEqual(self.sequences, result)

    def test_pattern_replace_static(self):
        result = transform.name_replace(self.sequences, '_REPLACE_',
                '_DONE_')
        result = list(result)
        expected = self.create_sequences()
        expected[1].id = 'test_DONE_2'
        self.assertEqual(self.sequences, result)

    def test_pattern_replace_case_insensitive(self):
        """
        Substitutions are case insensitive
        """
        result = transform.name_replace(self.sequences, '_replace_',
                '_DONE_')
        result = list(result)
        expected = self.create_sequences()
        expected[1].id = 'test_DONE_2'
        self.assertEqual(self.sequences, result)

    def test_pattern_replace_group(self):
        """
        Make sure capturing groups work
        """
        result = transform.name_replace(self.sequences, '_(repl)ace_',
                '_DONE-\\1_')
        result = list(result)
        expected = self.create_sequences()
        expected[1].id = 'test_DONE-repl_2'
        self.assertEqual(self.sequences, result)

class SqueezeTestCase(unittest.TestCase):

    def setUp(self):
        super(SqueezeTestCase, self).setUp()

        self.sequences = [
            seqrecord('sequence_1', 'AC-G--'),
            seqrecord('sequence_2', '-C-GT-'),
            seqrecord('sequence_3', '-T-AG-'),
        ]

    def test_gap_proportion(self):
        actual = transform.gap_proportion(self.sequences)
        self.assertEqual([2./3, 0.0, 1.0, 0.0, 1./3, 1.0], actual)

    def test_basic_squeeze(self):
        result = list(transform.squeeze(self.sequences, 1.0))

        self.assertEqual([4, 4, 4], [len(i) for i in result])
        self.assertEqual([i.id for i in self.sequences], [i.id for i in result])
        expected = [
            seqrecord('sequence_1', 'ACG-'),
            seqrecord('sequence_2', '-CGT'),
            seqrecord('sequence_3', '-TAG'),
        ]

        self.assertEqual([str(i.seq) for i in expected],
                [str(i.seq) for i in result])

    def test_squeeze_none(self):
        """
        Threshold of 0.001 - nothing should be squeezed.
        """
        result = list(transform.squeeze(self.sequences, 1.01))
        self.assertEqual([str(i.seq) for i in self.sequences],
                [str(i.seq) for i in result])


class SeqPatternTestCase(unittest.TestCase):

    def setUp(self):
        super(SeqPatternTestCase, self).setUp()

        self.sequences = [
            seqrecord('sequence_1', 'AC-G--'),
            seqrecord('sequence_2', '-C-GT-'),
            seqrecord('sequence_3', '-T-AG-'),
        ]

        self.tests = [('^$', []), ('.*', self.sequences),
                ('^AC', [self.sequences[0]])]

    def test_include(self):
        result = transform.seq_include(self.sequences, '^$')

        for regex, expected in self.tests:
            result = list(transform.seq_include(self.sequences, regex))
            self.assertEqual(expected, result)

    def test_exclude(self):
        result = transform.seq_include(self.sequences, '^$')

        for regex, expected_include in self.tests:
            expected = [i for i in self.sequences if i not in expected_include]
            result = list(transform.seq_exclude(self.sequences, regex))
            self.assertEqual(expected, result)

class HeadTestCase(unittest.TestCase):
    """
    Test for transform.head
    """

    def setUp(self):
        self.sequences = [seqrecord('sequence{0}'.format(i), 'A'*(i+1))
                          for i in xrange(100)]

    def test_zero(self):
        result = list(transform.head(self.sequences, 0))
        self.assertEqual([], result)

    def test_more_seqs_than_available(self):
        """
        Specifying more sequences than are in input records should return
        them all
        """
        result = list(transform.head(self.sequences, 10000))
        self.assertEqual(self.sequences, result)

    def test_values(self):
        """
        Try specifying some values.
        """
        for h in xrange(len(self.sequences) + 1):
            result = list(transform.head(self.sequences, h))
            self.assertEqual(h, len(result))
            self.assertEqual(self.sequences[:h], result)

class TailTestCase(unittest.TestCase):
    def setUp(self):
        self.records = [
            seqrecord('sequence_1', 'AC-G--'),
            seqrecord('sequence_2', '-C-GT-'),
            seqrecord('sequence_3', '-T-AG-'),
        ]

    def _do_test(self, size):
        actual = list(transform.tail(self.records, size))
        expected = self.records[-size:]
        self.assertEqual([e.id for e in expected], [a.id for a in actual])
        self.assertEqual([str(e.seq) for e in expected], [str(a.seq) for a in actual])

    def test_tail_1(self):
        self._do_test(1)

    def test_tail_2(self):
        self._do_test(2)

    def test_tail_3(self):
        self._do_test(3)

class IsolateRegionTestCase(unittest.TestCase):

    def setUp(self):
        self.sequences = [_alignment_record('--A--ACTGGACGTATTC-CCCC'),
                          _alignment_record('--AGCACTGGA---ATTC-CCCC')]

    def test_no_isolation(self):
        result = list(transform.isolate_region(self.sequences, 0,
            len(self.sequences[0])))

        self.assertEqual(self.sequences, result)

    def test_single_loc(self):
        start = 2
        end = 3
        result = list(transform.isolate_region(self.sequences, start, end))
        for seq in result:
            self.assertEqual('--A--------------------', str(seq.seq))

    def test_middle(self):
        expected = ['--A--ACTGGA------------', '--AGCACTGGA------------']
        start = 1
        end = 11

        actual = list(transform.isolate_region(self.sequences, start, end))
        actual = [str(s.seq) for s in actual]
        self.assertEqual(expected, actual)

    def test_invalid(self):
        self.assertRaises(ValueError, transform.isolate_region(
                self.sequences, 5, 5).next)
        self.assertRaises(ValueError, transform.isolate_region(
                self.sequences, 10, 5).next)

class MinUngapLengthTestCase(unittest.TestCase):

    def setUp(self):
        self.sequences = [_alignment_record('--AAC--'),
                          _alignment_record('AAAA...'),
                          _alignment_record('-------'),
                          _alignment_record('ACGRAGT')]

    def test_none_pass(self):
        result = list(transform.min_ungap_length_discard(self.sequences, 8))
        self.assertEqual([], result)

    def test_all_pass(self):
        result = list(transform.min_ungap_length_discard(self.sequences, 0))
        self.assertEqual(self.sequences, result)

    def test_partial(self):
        result = transform.min_ungap_length_discard(self.sequences, 4)
        self.assertEqual([self.sequences[1], self.sequences[3]], list(result))

class IncludeExcludeMixIn(object):

    def setUp(self):
        ids = """sequenceid1
sequenceid2
sequenceid4
"""
        self.handle = StringIO(ids)

        self.sequences = [SeqRecord(Seq("AAA"), id="sequenceid1"),
                SeqRecord(Seq("BBB"), id="sequenceid2"),
                SeqRecord(Seq("CCC"), id="sequenceid3"),
                SeqRecord(Seq("DDD"), id="sequenceid4", description='sequence id 4'),
                SeqRecord(Seq("EEE"), id="test", description='test sequence'), ]


class IncludeFromFileTestCase(IncludeExcludeMixIn, unittest.TestCase):

    def test_filter(self):
        expected = [self.sequences[0], self.sequences[1], self.sequences[3]]
        actual = list(transform.include_from_file(self.sequences, self.handle))
        self.assertEqual(3, len(actual))
        self.assertEqual(expected, actual)

class ExcludeFromFileTestCase(IncludeExcludeMixIn, unittest.TestCase):

    def test_filter(self):
        expected = [self.sequences[2], self.sequences[4]]
        actual = list(transform.exclude_from_file(self.sequences, self.handle))
        self.assertEqual(2, len(actual))
        self.assertEqual(expected, actual)

class NameIncludeTestCase(IncludeExcludeMixIn, unittest.TestCase):

    def test_filter_id(self):
        expected = self.sequences[:2]
        actual = list(transform.name_include(self.sequences, r'sequenceid[12]'))
        self.assertEqual(2, len(actual))
        self.assertEqual(expected, actual)

    def test_filter_description(self):
        expected = self.sequences[3:]
        actual = list(transform.name_include(self.sequences, r'sequence id 4|test seq'))
        self.assertEqual(2, len(actual))
        self.assertEqual(expected, actual)

class NameExcludeTestCase(IncludeExcludeMixIn, unittest.TestCase):

    def test_filter_id(self):
        expected = self.sequences[2:]
        actual = list(transform.name_exclude(self.sequences, r'sequenceid[12]'))
        self.assertEqual(3, len(actual))
        self.assertEqual(expected, actual)

    def test_filter_description(self):
        expected = self.sequences[:3]
        actual = list(transform.name_exclude(self.sequences, r'sequence id 4|test seq'))
        self.assertEqual(expected, actual)

class CutTestCase(unittest.TestCase):

    def setUp(self):
        self.sequences = [SeqRecord(Seq("ABC"), id="sequenceid1"),
                SeqRecord(Seq("BCD"), id="sequenceid2"),
                SeqRecord(Seq("DEF"), id="sequence id 4"),
                SeqRecord(Seq("EFG"), id="test sequence"), ]

    def test_no_sequences(self):
        actual = list(transform._cut_sequences(self.sequences, slice(0, 0)))
        for sequence in actual:
            self.assertEqual(0, len(sequence))

    def test_full_sequence(self):
        actual = list(transform._cut_sequences(self.sequences, slice(0, 3)))
        self.assertEqual(['ABC', 'BCD', 'DEF', 'EFG'], [str(s.seq) for s in
            actual])

    def test_cut_sequences(self):
        actual = list(transform._cut_sequences(self.sequences, slice(0, 2)))
        self.assertEqual(['AB', 'BC', 'DE', 'EF'], [str(s.seq) for s in
            actual])
        actual = list(transform._cut_sequences(self.sequences, slice(1, None)))
        self.assertEqual(['BC', 'CD', 'EF', 'FG'], [str(s.seq) for s in
            actual])

class CodonWarningTableTestCase(unittest.TestCase):

    def warn(self, *args, **kwargs):
        self.warnings.append((args, kwargs))

    def setUp(self):
        self.warnings = []
        self.warning_dict = transform.CodonWarningTable({'UUU': 'F'})
        self.old_warn = transform.logging.warn
        transform.logging.warn = self.warn

    def tearDown(self):
        transform.logging.warn = self.old_warn

    def test_nowarn(self):
        actual = self.warning_dict['UUU']
        self.assertEqual('F', actual)
        self.assertEqual([], self.warnings)

    def test_warn(self):
        codon = 'UU-'
        actual = self.warning_dict[codon]
        self.assertEqual('X', actual)
        self.assertEqual([(("Unknown Codon: %s", codon), {})], self.warnings)

class TranslateTestCase(unittest.TestCase):

    def test_dna_protein_nogap(self):
        sequences = [seqrecord('A', 'TTTTTATAA')]
        expected = ['FL*']
        actual = transform.translate(sequences, 'dna2protein')
        self.assertEqual(expected, [str(i.seq) for i in actual])

    def test_dna_protein_nogap_stop(self):
        sequences = [seqrecord('A', 'TTTTTATAA')]
        expected = ['FL']
        actual = transform.translate(sequences, 'dna2proteinstop')
        self.assertEqual(expected, [str(i.seq) for i in actual])

    def test_dna_protein_gap(self):
        sequences = [seqrecord('A', 'TTTTT-TAA')]
        expected = ['FX*']
        actual = transform.translate(sequences, 'dna2protein')
        self.assertEqual(expected, [str(i.seq) for i in actual])

    def test_dna_protein_gap_stop(self):
        sequences = [seqrecord('A', '---TTATAA')]
        expected = ['-L']
        actual = transform.translate(sequences, 'dna2proteinstop')
        self.assertEqual(expected, [str(i.seq) for i in actual])

class UngapSequencesTestCase(unittest.TestCase):

    def test_dot_gap(self):
        sequences = [SeqRecord(Seq("AAA"), id="s1"),
                SeqRecord(Seq("A.G"), id="s2"),
                SeqRecord(Seq(".A."), id="s3"),]

        ungapped = list(transform.ungap_sequences(sequences))
        self.assertEqual(["AAA", "AG", "A"], [str(s.seq) for s in ungapped])

    def test_dash_gap(self):
        sequences = [SeqRecord(Seq("AAA"), id="s1"),
                SeqRecord(Seq("A-G"), id="s2"),
                SeqRecord(Seq("-A-"), id="s3"),]

        ungapped = list(transform.ungap_sequences(sequences))
        self.assertEqual(["AAA", "AG", "A"], [str(s.seq) for s in ungapped])

# Name Modification functions
class IdModifyMixin(object):
    """
    Mixin to ease testing name prefix and suffix
    """

    def setUp(self):
        self.input_fp = StringIO(self.initial_fasta)
        self.output_fp = StringIO()

    def test_modify(self):
        records = SeqIO.parse(self.input_fp, 'fasta')
        records = self.__class__.modify_fn(records)
        SeqIO.write(records, self.output_fp, 'fasta')
        self.assertEqual(self.target_fasta, self.output_fp.getvalue().strip())

class NamePrefixTestCase(IdModifyMixin, unittest.TestCase):
    initial_fasta = """>seq1
ACGT
>gi|260674|gb|S52561.1| {long terminal repeat} [human immunodeficiency virus type]
ACGT"""
    target_fasta = """>pre.seq1
ACGT
>pre.gi|260674|gb|S52561.1| {long terminal repeat} [human immunodeficiency virus type]
ACGT"""
    modify_fn = functools.partial(transform.name_insert_prefix, prefix="pre.")

class NameSuffixTestCase(IdModifyMixin, unittest.TestCase):
    initial_fasta = """>seq1
ACGT
>gi|260674|gb|S52561.1| {long terminal repeat} [human immunodeficiency virus type]
ACGT"""
    target_fasta = """>seq1.post
ACGT
>gi|260674|gb|S52561.1|.post {long terminal repeat} [human immunodeficiency virus type]
ACGT"""
    modify_fn = functools.partial(transform.name_append_suffix, suffix=".post")


class MultiCutTestCase(unittest.TestCase):
    def setUp(self):
        self.inputs = [seqrecord("Sequence 1", "ACGT--TCAGA")]

    def test_multicut(self):
        actual = list(transform.multi_cut_sequences(self.inputs,
            [slice(None, 2), slice(8, None)]))
        self.assertEqual(['ACAGA'], [str(s.seq) for s in actual])

class MultiMaskSequences(unittest.TestCase):

    def setUp(self):
        self.sequences = [SeqRecord(Seq("AAA"), id="sequenceid1"),
                SeqRecord(Seq("BBB"), id="sequenceid2"),
                SeqRecord(Seq("DDDD"), id="sequence id 4"),
                SeqRecord(Seq("EEE"), id="test sequence"), ]

    def test_mask_whole(self):
        masks = [slice(0, 200)]
        actual = list(transform.multi_mask_sequences(self.sequences, masks))
        self.assertEqual(len(self.sequences), len(actual))
        for e, a in zip(self.sequences, actual):
            self.assertEqual(e.id, a.id)
            self.assertEqual('-'*len(e), str(a.seq))

    def test_mask(self):
        masks = [slice(1, 2)]
        actual = list(transform.multi_mask_sequences(self.sequences, masks))
        self.assertEqual(len(self.sequences), len(actual))
        self.assertEqual(['A-A', 'B-B', 'D-DD', 'E-E'],
                [str(a.seq) for a in actual])

class RecordBufferTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = [SeqRecord(Seq("AAA"), id="s1"),
                SeqRecord(Seq("A-G"), id="s2"),
                SeqRecord(Seq("-A-"), id="s3"),]
        self.seq_iter = iter(self.sequences)

    def _compare(self, records):
        self.assertEqual(len(self.sequences), len(records))

        for e, a in zip(self.sequences, records):
            self.assertEqual(e.id, a.id)
            self.assertEqual(e.description, a.description)
            self.assertEqual(str(e.seq), str(a.seq))

    def test_single_pass(self):
        with transform._record_buffer(self.seq_iter) as iter_f:
            records = list(iter_f())
            self._compare(records)

    def test_multi_pass(self):
        with transform._record_buffer(self.seq_iter) as iter_f:
            records = list(iter_f())
            self._compare(records)

            records = list(iter_f())
            self._compare(records)

class DropColumnsTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = [SeqRecord(Seq("AAA"), id="s1"),
                SeqRecord(Seq("A-G"), id="s2"),
                SeqRecord(Seq("-A-"), id="s3"),]

    def test_basic(self):
        r = list(transform.drop_columns(self.sequences, [slice(1, None)]))
        self.assertEqual([i.id for i in self.sequences],
                [i.id for i in r])
        self.assertEqual(['A', 'A', '-'], [str(i.seq) for i in r])

    def test_multi(self):
        r = list(transform.drop_columns(self.sequences, [slice(0, 1), slice(2, None)]))
        self.assertEqual([i.id for i in self.sequences],
                [i.id for i in r])
        self.assertEqual(['A', '-', 'A'], [str(i.seq) for i in r])

class DashesCleanupTestCase(unittest.TestCase):
    def setUp(self):
        self.sequences = [SeqRecord(Seq("A~-.?~GT"), id="s1"),
                SeqRecord(Seq("A-GGGG?-"), id="s2"),
                SeqRecord(Seq("-A-:ACA-"), id="s3"),
                SeqRecord(Seq("ACTGGTCA"), id="s4"),]

    def test_basic(self):
        actual = list(transform.dashes_cleanup(self.sequences))
        actual = [(i.id, str(i.seq)) for i in actual]
        self.assertEqual(
                [('s1', 'A-----GT'),
                 ('s2', 'A-GGGG--'),
                 ('s3', '-A--ACA-'),
                 ('s4', 'ACTGGTCA')], actual)
