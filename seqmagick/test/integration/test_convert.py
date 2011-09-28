from cStringIO import StringIO
import os
import os.path
import shlex
import shutil
import sys
import unittest
import tempfile

from seqmagick.scripts import cli


d = os.path.dirname(__file__)
data_dir = os.path.join(d, "data")

def p(*args):
    return os.path.join(data_dir, *args)

class CommandLineTestMixIn(object):
    in_suffix = ''
    out_suffix = ''

    def setUp(self):
        self.input_file = tempfile.NamedTemporaryFile(suffix=self.in_suffix)
        with open(self.input_path) as fp:
            shutil.copyfileobj(fp, self.input_file)
        self.input_file.flush()
        with tempfile.NamedTemporaryFile(suffix=self.out_suffix) as tf:
            self.output_file = tf.name

    def test_run(self):
        command = self.command.format(input=self.input_file.name,
                output=self.output_file)
        cli.main(shlex.split(command))

        with open(self.output_file) as fp:
            actual = fp.read()
        with open(self.expected_path) as fp:
            expected = fp.read()
        self.assertEqual(expected, actual)

    def tearDown(self):
        self.input_file.close()
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

class TestBasicConvert(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.phy'
    input_path = p('input2.fasta')
    expected_path = p('output2.phy')
    command = 'convert {input} {output}'

class TestConvertUngapCut(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.fasta'
    input_path = p('input2.fasta')
    expected_path = p('output2_ungap_cut.fasta')
    command = 'convert --ungap --cut 1:3 --tail 2 {input} {output}'

class TestConvertToStdOut(unittest.TestCase):

    def setUp(self):
        self.out = StringIO()
        self.actual_stdout = sys.stdout
        sys.stdout = self.out

    def tearDown(self):
        sys.stdout = self.actual_stdout

    def test_convert(self):
        in_path = p('input2.fasta')
        cli.main(['convert', in_path, '-', '--output-format', 'fasta'])
        actual = self.out.getvalue()
        with open(in_path) as fp:
            expected = fp.read()
        self.assertEqual(expected, actual)
