import os.path
import unittest
import tempfile

from seqmagick.scripts import cli

d = os.path.dirname(__file__)
data_dir = os.path.join(d, "data")

class TestInfo(unittest.TestCase):

    def setUp(self):
        self.infile = tempfile.NamedTemporaryFile()
        self.tempfile = tempfile.NamedTemporaryFile()

    def tearDown(self):
        self.infile.close()
        self.tempfile.close()

    def test_simple(self):
        seq_file = os.path.join(data_dir, 'input1.fasta')
        args = ['info', seq_file,
                '--out-file', self.tempfile.name]

        cli.main(args)
        self.assertEquals("""name\talignment\tmin_len\tmax_len\tavg_len\tnum_seqs
{0}\tTRUE\t4\t4\t4.00\t3
""".format(seq_file), self.tempfile.read())

