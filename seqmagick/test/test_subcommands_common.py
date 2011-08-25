import argparse
import unittest

from seqmagick.subcommands import common

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
        self.assertEqual([0, 1], f(range(10)))

    def test_no_arg(self):
        def head(records):
            return records[:2]

        a_cls = common.partial_append_action(head)
        a = a_cls([], 'functions')

        a(None, self.namespace, None)
        f = self.namespace.functions[0]
        self.assertEqual([0, 1], f(range(10)))

    def test_multi_arg(self):
        def fake_slice(records, i, j):
            return records[i:j]

        a_cls = common.partial_append_action(fake_slice, ['i', 'j'])
        a = a_cls([], 'functions')

        a(None, self.namespace, [0, 2])
        f = self.namespace.functions[0]
        self.assertEqual([0, 1], f(range(10)))

class PositiveValueTestCase(unittest.TestCase):

    def test_negative(self):
        self.assertRaises(argparse.ArgumentTypeError,
                common.positive_value(int), '-1')

    def test_positive(self):
        self.assertEqual(1, common.positive_value(int)('1'))

    def test_zero(self):
        self.assertEqual(0, common.positive_value(int)('0'))

class CutRangeTestCase(unittest.TestCase):
    def test_negative(self):
        self.assertRaises(argparse.ArgumentTypeError,
                common.cut_range, '0:5')

    def test_out_of_order(self):
        self.assertRaises(argparse.ArgumentTypeError,
                common.cut_range, '10:5')

    def test_start(self):
        actual = common.cut_range('5:10')
        self.assertEqual(4, actual.start)
        self.assertEqual(10, actual.stop)
