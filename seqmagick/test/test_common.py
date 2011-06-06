import argparse
import unittest

from seqmagick.subcommands import common


class PositiveValueTestCase(unittest.TestCase):

    def test_negative(self):
        self.assertRaises(argparse.ArgumentTypeError,
                common.positive_value(int), '-1')

    def test_positive(self):
        self.assertEquals(1, common.positive_value(int)('1'))

    def test_zero(self):
        self.assertEquals(0, common.positive_value(int)('0'))
