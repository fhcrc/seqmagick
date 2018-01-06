import argparse
import os
import os.path
import sys
import unittest
import tempfile

from seqmagick.subcommands import common

d = os.path.dirname(__file__)
data_dir = os.path.join(d, "integration", "data")


def p(*args):
    return os.path.join(data_dir, *args)


class PartialAppendTestCase(unittest.TestCase):

    def setUp(self):
        self.namespace = argparse.Namespace()

    def test_single_arg(self):
        def head(records, n):
            return records[:n]

        a_cls = common.partial_append_action(head, 'n')
        a = a_cls([], 'functions')

        a(None, self.namespace, 2)

        f = self.namespace.functions[0]
        self.assertEqual([0, 1], f(list(range(10))))

    def test_no_arg(self):
        def head(records):
            return records[:2]

        a_cls = common.partial_append_action(head)
        a = a_cls([], 'functions')

        a(None, self.namespace, None)
        f = self.namespace.functions[0]
        self.assertEqual([0, 1], f(list(range(10))))

    def test_multi_arg(self):
        def fake_slice(records, i, j):
            return records[i:j]

        a_cls = common.partial_append_action(fake_slice, ['i', 'j'])
        a = a_cls([], 'functions')

        a(None, self.namespace, [0, 2])
        f = self.namespace.functions[0]
        self.assertEqual([0, 1], f(list(range(10))))

class PositiveValueTestCase(unittest.TestCase):

    def test_negative(self):
        self.assertRaises(argparse.ArgumentTypeError,
                common.positive_value(int), '-1')

    def test_positive(self):
        self.assertEqual(1, common.positive_value(int)('1'))

    def test_zero(self):
        self.assertEqual(0, common.positive_value(int)('0'))

class CutRangeTestCase(unittest.TestCase):
    def test_out_of_order(self):
        self.assertRaises(argparse.ArgumentTypeError,
                common.cut_range, '10:5')

    def test_start(self):
        actual = common.cut_range('5:10')
        self.assertEqual(4, actual.start)
        self.assertEqual(10, actual.stop)

    def test_negative(self):
        actual = common.cut_range('-500:')
        self.assertEqual(-500, actual.start)
        self.assertIsNone(actual.stop)
        actual = common.cut_range('-500:-203')
        self.assertEqual(-500, actual.start)
        self.assertEqual(-203, actual.stop)

    def test_no_start(self):
        actual = common.cut_range(':10')
        self.assertEqual(None, actual.start)
        self.assertEqual(10, actual.stop)

    def test_no_end(self):
        actual = common.cut_range('5:')
        self.assertEqual(4, actual.start)

class SequenceSlicesTestCase(unittest.TestCase):
    def test_single(self):
        actual = common.sequence_slices(':10')
        self.assertEqual([slice(None, 10)], actual)

    def test_multiple(self):
        actual = common.sequence_slices('1:10,3:20')
        self.assertEqual([slice(0, 10), slice(2, 20)], actual)

class AtomicWriteTestCase(unittest.TestCase):

    initial_content = "Initial Content"
    new_content = "New Content"

    def setUp(self):
        with tempfile.NamedTemporaryFile('wt', delete=False) as tf:
            tf.write(self.initial_content)
            self.input_file = tf.name

    def test_exception_leaves_unchanged(self):
        try:
            with common.atomic_write(self.input_file) as tf:
                raise IOError()
        except IOError:
            with open(self.input_file) as fp:
                self.assertEqual(self.initial_content, fp.read())

            # Ensure deleted
            self.assertFalse(os.path.exists(tf.name))

    def test_write(self):
        with common.atomic_write(self.input_file) as fp:
            self.assertNotEqual(self.input_file, fp.name)
            fp.write(self.new_content)

        self.assertFalse(os.path.exists(fp.name))

        with open(self.input_file) as fp:
            self.assertEqual(self.new_content, fp.read())

    def tearDown(self):
        os.remove(self.input_file)


class ApplyUmaskTestCase(unittest.TestCase):

    def setUp(self):
        # Set umask
        self.orig_umask = common.get_umask()

    def tearDown(self):
        os.umask(self.orig_umask)

    def test_provided_umask(self):
        self.assertEqual('0o770', oct(common.apply_umask(0o777, 0o07)))
        self.assertEqual('0o660', oct(common.apply_umask(0o666, 0o07)))
        self.assertEqual('0o644', oct(common.apply_umask(0o666, 0o22)))

    def test_user_umask(self):
        os.umask(0o07)
        self.assertEqual('0o770', oct(common.apply_umask(0o777)))
        self.assertEqual('0o660', oct(common.apply_umask(0o666)))


class FileTypeTestCase(unittest.TestCase):

    def setUp(self):
        # used in methods test_read_*()
        self.testfile = 'input2.fasta'
        with open(p(self.testfile)) as f:
            self.expected = f.read()

    def test_stdin(self):
        self.assertIs(sys.stdin, common.FileType('r')('-'))

    def test_stdout(self):
        self.assertIs(sys.stdout, common.FileType('w')('-'))

    def test_read(self):
        with tempfile.NamedTemporaryFile('w+t') as tf:
            tf.write('TEST')
            tf.flush()
            with common.FileType('r')(tf.name) as fp:
                self.assertEqual(tf.name, fp.name)
                self.assertEqual('TEST', fp.read())

    def test_write(self):
        with tempfile.NamedTemporaryFile('w+t') as tf:
            with common.FileType('w')(tf.name) as fp:
                fp.write('TEST')
                fp.flush()
                self.assertEqual(tf.name, fp.name)
                self.assertEqual('TEST', tf.read())

    def test_read_text(self):
        with common.FileType('rt')(p(self.testfile)) as fp:
            self.assertEqual(fp.read(), self.expected)

    def test_read_gz(self):
        with common.FileType('rt')(p(self.testfile + '.gz')) as fp:
            self.assertEqual(fp.read(), self.expected)

    @unittest.skipIf(sys.version_info.major == 3, 'bzip2 not supported')
    def test_read_bz2(self):
        with common.FileType('rt')(p(self.testfile + '.bz2')) as fp:
            self.assertEqual(fp.read(), self.expected)
