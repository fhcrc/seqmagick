import unittest
import tempfile

from seqmagick.scripts import cli

from seqmagick.test.integration import data_path

class ExtractIdsMixin(object):
    expected = """test1
test2
test3
"""
    expected_desc = """test1 test sequence 1
test2 test sequence 2
test3 sequence 3
"""

    def setUp(self):
        self.tempfile = tempfile.NamedTemporaryFile()

    def tearDown(self):
        self.tempfile.close()

    def test_ids(self):
        args = ['extract-ids', self.seq_file, '-o', self.tempfile.name]
        cli.main(args)
        self.assertEquals(self.expected, self.tempfile.read())

    def test_descriptions(self):
        args = ['extract-ids', self.seq_file, '-o', self.tempfile.name, '-d']
        cli.main(args)
        self.assertEquals(self.expected_desc, self.tempfile.read())


class SimpleExtractIdsTestCase(ExtractIdsMixin, unittest.TestCase):
    seq_file = data_path('input2.fasta')

class Bz2ExtractIdsTestCase(ExtractIdsMixin, unittest.TestCase):
    seq_file = data_path('input2.fasta.bz2')

class GzipExtractIdsTestCase(ExtractIdsMixin, unittest.TestCase):
    seq_file = data_path('input2.fasta.gz')
