
import os
import os.path
import shlex
import shutil
import tempfile

from seqmagick.scripts import cli

from . import test_convert

d = os.path.dirname(__file__)
data_dir = os.path.join(d, "data")

def p(*args):
    return os.path.join(data_dir, *args)

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

        with open(self.input_file) as fp:
            actual = fp.read()
        with open(self.expected_path) as fp:
            expected = fp.read()
        self.assertEqual(expected, actual)

    def tearDown(self):
        os.remove(self.input_file)

class TestMogrifyUngapCut(CommandLineTestMixIn, test_convert.TestConvertUngapCut):
    command = 'mogrify --ungap --cut 1:3 --tail 2 {input}'
