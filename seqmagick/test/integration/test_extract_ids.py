import os.path
import unittest
import tempfile

from seqmagick.scripts import cli


d = os.path.dirname(__file__)
data_dir = os.path.join(d, "data")

class TestExtractIds(unittest.TestCase):

    def setUp(self):
        self.tempfile = tempfile.NamedTemporaryFile()

    def tearDown(self):
        self.tempfile.close()

    def test_simple(self):
        args = ['extract-ids', os.path.join(data_dir, 'input1.fasta'),
                '-o', self.tempfile.name]
        cli.main(args)
        self.assertEquals("""test1
test2
test3
""", self.tempfile.read())

    def test_description(self):
        args = ['extract-ids', os.path.join(data_dir, 'input1.fasta'),
                '-o', self.tempfile.name, '-d']
        cli.main(args)
        self.assertEquals("""test1 test sequence 1
test2 test sequence 2
test3 sequence 3
""", self.tempfile.read())
