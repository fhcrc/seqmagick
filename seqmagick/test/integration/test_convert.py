from cStringIO import StringIO
import os
import os.path
import logging
import shlex
import shutil
import sys
import unittest
import tempfile

from seqmagick.subcommands.common import FileType
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

        with FileType('r')(self.output_file) as fp:
            actual = fp.read()
        with FileType('r')(self.expected_path) as fp:
            expected = fp.read()
        self.assertEqual(expected, actual)

    def tearDown(self):
        self.input_file.close()
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

class BasicConvertTestCase(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.phy'
    input_path = p('input2.fasta')
    expected_path = p('output2.phy')
    command = 'convert {input} {output}'

class BzipInputConvertTestCase(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta.bz2'
    out_suffix = '.phy'
    input_path = p('input2.fasta.bz2')
    expected_path = p('output2.phy')
    command = 'convert {input} {output}'

class BzipOutputConvertTestCase(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.phy.bz2'
    input_path = p('input2.fasta')
    expected_path = p('output2.phy')
    command = 'convert {input} {output}'

class GzipInputConvertTestCase(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta.gz'
    out_suffix = '.phy'
    input_path = p('input2.fasta.gz')
    expected_path = p('output2.phy')
    command = 'convert {input} {output}'

class GzipOutputConvertTestCase(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.phy.gz'
    input_path = p('input2.fasta')
    expected_path = p('output2.phy')
    command = 'convert {input} {output}'

class ConvertToNexusTestCase(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    input_path = p('input2.fasta')
    expected_path = p('output2.nex')
    command = 'convert {input} {output} --output-format nexus --alphabet dna-ambiguous'

class ConvertUngapCutTestCase(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.fasta'
    input_path = p('input2.fasta')
    expected_path = p('output2_ungap_cut.fasta')
    command = 'convert --ungap --cut 1:3 --tail 2 {input} {output}'

class ConvertToStdOutTestCase(unittest.TestCase):

    def setUp(self):
        self.out = StringIO()
        self.err = StringIO()
        self.actual_stdout = sys.stdout
        self.actual_stderr = sys.stderr
        sys.stdout = self.out
        sys.stderr = self.err

    def tearDown(self):
        sys.stdout = self.actual_stdout
        sys.stderr = self.actual_stderr

    def test_convert(self):
        in_path = p('input2.fasta')
        cli.main(['convert', in_path, '-', '--output-format', 'fasta'])
        actual = self.out.getvalue()
        with open(in_path) as fp:
            expected = fp.read()
        self.assertEqual(expected, actual)

    def test_convert_noformat(self):
        in_path = p('input2.fasta')
        cli.main(['convert', in_path, '-'])
        actual = self.out.getvalue()
        with open(in_path) as fp:
            expected = fp.read()
        self.assertEqual(expected, actual)

class TestCutRelative(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.fasta'
    input_path = p('input3.fasta')
    expected_path = p('output3.fasta')
    command = 'convert --cut 2:3 --relative-to HXB2 {input} {output}'

    def test_unknown_seq(self):
        args = ['convert', '--cut', '2:3', '--relative-to', 'OTHER',
                self.input_path, '-', '--output-format', 'fasta']
        self.assertRaises(ValueError, cli.main, args)

class TestTranslateAmbiguous(CommandLineTestMixIn, unittest.TestCase):
    in_suffix = '.fasta'
    out_suffix = '.fasta'
    input_path = p('input4_ambig.fasta')
    expected_path = p('output4.fasta')
    command = 'convert --translate dna2protein {input} {output}'

    def setUp(self):
        super(TestTranslateAmbiguous, self).setUp()
        self.orig_level = logging.getLogger(None).level
        logging.getLogger(None).setLevel(logging.FATAL)

    def tearDown(self):
        super(TestTranslateAmbiguous, self).tearDown()
        logging.getLogger(None).setLevel(self.orig_level)

class TestStdin(TestTranslateAmbiguous, unittest.TestCase):
    command = 'convert --translate dna2protein - {output}'

    def setUp(self):
        super(TestStdin, self).setUp()
        self.orig_stdin = sys.stdin
        sys.stdin = open(p('input4_ambig.fasta'), 'r')

    def tearDown(self):
        super(TestStdin, self).tearDown()
        sys.stdin.close()
        sys.stdin = self.orig_stdin

class TestConvertFromStdin(TestTranslateAmbiguous, unittest.TestCase):
    command = 'convert --translate dna2protein - {output}'

    def setUp(self):
        super(TestConvertFromStdin, self).setUp()
        self.orig_stdin = sys.stdin
        sys.stdin = open(p('input4_ambig.fasta'), 'r')

    def tearDown(self):
        super(TestConvertFromStdin, self).tearDown()
        sys.stdin.close()
        sys.stdin = self.orig_stdin
