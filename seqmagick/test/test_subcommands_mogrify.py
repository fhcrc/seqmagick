import os
import os.path
import tempfile
import unittest

from seqmagick.subcommands import mogrify

class AtomicWriteTestCase(unittest.TestCase):

    initial_content = "Initial Content"
    new_content = "New Content"

    def setUp(self):
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            tf.write(self.initial_content)
            self.input_file = tf.name

    def test_exception_leaves_unchanged(self):
        try:
            with mogrify.atomic_write(self.input_file) as tf:
                raise IOError()
        except IOError:
            with open(self.input_file) as fp:
                self.assertEqual(self.initial_content, fp.read())

            # Ensure deleted
            self.assertFalse(os.path.exists(tf.name))

    def test_write(self):
        with mogrify.atomic_write(self.input_file) as fp:
            self.assertNotEqual(self.input_file, fp.name)
            fp.write(self.new_content)

        self.assertFalse(os.path.exists(fp.name))

        with open(self.input_file) as fp:
            self.assertEqual(self.new_content, fp.read())

    def tearDown(self):
        os.remove(self.input_file)
