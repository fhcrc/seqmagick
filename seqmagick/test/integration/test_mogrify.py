
import os
import os.path
import shlex
import shutil
import tempfile

from seqmagick.scripts import cli
from seqmagick.subcommands.common import FileType
from seqmagick.test.integration import data_path

from . import test_convert

class CommandLineTestMixIn(object):
    def setUp(self):
        with tempfile.NamedTemporaryFile(
                suffix=os.path.basename(self.input_path),
                delete=False) as tf:
            self.input_file = tf.name
            with open(self.input_path) as fp:
                shutil.copyfileobj(fp, tf)

    def test_run(self):
        command = self.command.format(input=self.input_file)
        try:
            cli.main(shlex.split(command))
        except SystemExit, e:
            self.fail(e)

        with FileType('r')(self.input_file) as fp:
            actual = fp.read()
        with FileType('r')(self.expected_path) as fp:
            expected = fp.read()
        self.assertEqual(expected, actual)

    def tearDown(self):
        os.remove(self.input_file)

class MogrifyUngapCutTestCase(CommandLineTestMixIn, test_convert.ConvertUngapCutTestCase):
    command = 'mogrify --ungap --cut 1:3 --tail 2 {input}'

class MogrifyBzipInputTestCase(CommandLineTestMixIn, test_convert.BzipInputConvertTestCase):
    command = 'mogrify {input}'
    expected_path = data_path('output2.fasta')
    out_suffix='fasta.bz2'

class MogrifyGzipInputTestCase(CommandLineTestMixIn, test_convert.GzipInputConvertTestCase):
    command = 'mogrify {input}'
    expected_path = data_path('output2.fasta')
    out_suffix='fasta.gz'
