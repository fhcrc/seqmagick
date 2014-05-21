"""
Tests for primer trim
"""
import unittest

from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from seqmagick.subcommands.protparam import ProtParamCalculator

class ProtParamTestCase(unittest.TestCase):

    def setUp(self):
        self.sequences = [
            SeqRecord(Seq('MGHFTEEDKATITSLWGKVNVEDAGGETLGRLLVVYPWTQRFF'
                          'DSFGNLSSASAIMGNPKVKAHGKKVLTSLGDAIKHLDDLKGTF'
                          'AQLSELHCDKLHVDPENFKLLGNVLVTVLAIHFGKEFTPEVQA'
                          'SWQKMVTAVASALSSRYH'), id='1'),
            SeqRecord(Seq('MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTY'
                          'FPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSD'
                          'LHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFL'
                          'ASVSTVLTSKYR'), id='2'),
            SeqRecord(Seq('MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFF'
                          'ESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTF'
                          'ATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQA'
                          'AYQKVVAGVANALAHKYH'), id='3'),
        ]

        self.masses = [('1', 16139.86), ('2', 15256.89), ('3', 15997.81)]
        self.pis = [('1', 6.644836), ('2', 8.716858), ('3', 6.744568)]


        self.bad_sequence = [
            SeqRecord(Seq('XXHFTEEDKATITS'), id='99'),
        ]
        self.bad_seq_mass = 1567.36
        self.bad_seq_pi = 4.651184

    def test_molecular_weight_unsorted(self):
        self.instance = ProtParamCalculator(self.sequences)
        stats = self.instance.stats
        expected = self.masses
        for i, s in enumerate(stats):
            self.assertEqual(stats[i][0].id, expected[i][0])
            self.assertAlmostEqual(stats[i][1], expected[i][1], 2)

    def test_isoelectric_point_unsorted(self):
        self.instance = ProtParamCalculator(self.sequences)
        stats = self.instance.stats
        expected = self.pis
        for i, s in enumerate(stats):
            self.assertEqual(stats[i][0].id, expected[i][0])
            self.assertAlmostEqual(stats[i][2], expected[i][1], 6)

    def test_molecular_weight_sorted_desc(self):
        self.instance = ProtParamCalculator(self.sequences, sort_on='mass',
                                            sort_ascending=False)
        stats = self.instance.stats
        # sort expected masses descending
        expected = sorted(self.masses, key=lambda item: item[1], reverse=True)
        for i, s in enumerate(stats):
            self.assertEqual(stats[i][0].id, expected[i][0])
            self.assertAlmostEqual(stats[i][1], expected[i][1], 2)

    def test_isoelectric_point_sorted_asc(self):
        #self.instance.sort_on = 'pi'
        #self.instance.sort_ascending = True
        #stats = self.instance.calculate(self.sequences)
        self.instance = ProtParamCalculator(self.sequences, sort_on='pi',
                                            sort_ascending=True)
        stats = self.instance.stats
        # sort expected pis ascending
        expected = sorted(self.pis, key=lambda item: item[1])
        for i, s in enumerate(stats):
            self.assertEqual(stats[i][0].id, expected[i][0])
            self.assertAlmostEqual(stats[i][2], expected[i][1], 6)

    def test_recalc_after_changing_parameters(self):
        self.instance = ProtParamCalculator(self.sequences, sort_on='mass',
                                            sort_ascending=False)
        self.instance.sort_on = 'pi'
        self.instance.sort_ascending = True
        stats = self.instance.calculate(self.sequences)
        # sort expected pis ascending
        expected = sorted(self.pis, key=lambda item: item[1])
        for i, s in enumerate(stats):
            self.assertEqual(stats[i][0].id, expected[i][0])
            self.assertAlmostEqual(stats[i][2], expected[i][1], 6)

    def test_catch_unknown_residue(self):
        with self.assertRaises(ValueError) as error:
            self.instance = ProtParamCalculator(self.bad_sequence)

    def test_allow_unknown_residues(self):
        self.instance = ProtParamCalculator(self.bad_sequence,
                                            allow_unknown_residues=True)
        stats = self.instance.stats
        self.assertAlmostEqual(stats[0][1], self.bad_seq_mass, 2)
        self.assertAlmostEqual(stats[0][2], self.bad_seq_pi, 6)