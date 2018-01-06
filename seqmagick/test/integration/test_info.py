import sys
import unittest
import tempfile

from seqmagick.scripts import cli

from seqmagick.test.integration import data_path


class InfoMixin(object):
    expected = """name\talignment\tmin_len\tmax_len\tavg_len\tnum_seqs
{0}\tTRUE\t5\t5\t5.00\t3
"""
    threads = 1

    def setUp(self):
        self.infile = tempfile.NamedTemporaryFile()
        self.tempfile = tempfile.NamedTemporaryFile('w+t')

    def tearDown(self):
        self.infile.close()
        self.tempfile.close()

    def test_info(self):
        args = ['info', self.seq_file,
                '--out-file', self.tempfile.name,
                '--threads', str(self.threads)]

        cli.main(args)
        self.assertEqual(self.expected.format(self.seq_file), self.tempfile.read())


class SimpleInfoTestCase(InfoMixin, unittest.TestCase):
    seq_file = data_path('input2.fasta')


class MultithreadedInfoTestCase(InfoMixin, unittest.TestCase):
    seq_file = data_path('input2.fasta')
    threads = 2


class SimpleGzipInfoTestCase(InfoMixin, unittest.TestCase):
    seq_file = data_path('input2.fasta.gz')


@unittest.skipIf(sys.version_info.major == 3, 'bzip2 not supported')
class SimpleBzip2InfoTestCase(InfoMixin, unittest.TestCase):
    seq_file = data_path('input2.fasta.bz2')
