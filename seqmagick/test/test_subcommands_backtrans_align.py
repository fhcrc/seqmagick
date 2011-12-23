import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable

from seqmagick.subcommands import backtrans_align

class BatchTestCase(unittest.TestCase):
    def test_no_input(self):
        i = []
        b = backtrans_align.batch(i, 1)
        self.assertRaises(StopIteration, next, b)

    def test_singletons(self):
        i = range(3)
        b = backtrans_align.batch(i, 1)
        self.assertEquals([[0], [1], [2]], list(b))

    def test_doubles(self):
        i = range(6)
        b = backtrans_align.batch(i, 2)
        self.assertEquals([[0, 1], [2, 3], [4, 5]], list(b))

    def test_partial(self):
        i = range(5)
        b = backtrans_align.batch(i, 2)
        self.assertEquals([[0, 1], [2, 3], [4]], list(b))


class AlignmentMapperTestCase(unittest.TestCase):
    def setUp(self):
        self.instance = backtrans_align.AlignmentMapper(CodonTable.unambiguous_dna_by_name['Standard'])

    def test_validate_valid(self):
        nucl = 'TTTAAG'
        prot = 'FK'

        self.assertTrue(self.instance._validate_translation(prot, nucl))

    def test_validate_invalid(self):
        nucl = 'AAGTTT'
        prot = 'KK'
        self.assertRaisesRegexp(ValueError, r'Codon TTT translates to F, not K',
                self.instance._validate_translation, prot, nucl)

    def test_map_alignment(self):
        nucl = [SeqRecord(Seq('AAGTTT'), id='1'), # KF
                SeqRecord(Seq('AAGGTCTTC'), id='2'), # KVF
                SeqRecord(Seq('GGGGTTTTT'), id='3')] # GVF
        prot = [SeqRecord(Seq('-K-F'), id='1'),
                SeqRecord(Seq('-KVF'), id='2'),
                SeqRecord(Seq('G-VF'), id='3')]

        result = self.instance.map_all(prot, nucl)
        result = [(s.id, str(s.seq)) for s in result]
        self.assertEqual([('1', '---AAG---TTT'),
                          ('2', '---AAGGTCTTC'),
                          ('3', 'GGG---GTTTTT')], result)

        pass
