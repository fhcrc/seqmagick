import argparse
import unittest

from seqmagick.subcommands import common


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
