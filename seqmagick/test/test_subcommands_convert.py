"""
Tests for seqmagick.subcommands.convert - mostly integration with
seqmagick.transform
"""
import argparse
import os
import tempfile
import unittest

from seqmagick.subcommands import convert
from seqmagick import transform

# Test populating the transform
class PopulateTransformsMixIn(object):
    """
    Tests that transforms list is populated
    """
    def setUp(self):
        self.parser = convert.build_parser(argparse.ArgumentParser())
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            self.infile = tf.name
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            self.outfile = tf.name

    def tearDown(self):
        os.remove(self.infile)
        os.remove(self.outfile)

    def test_parse(self):
        arguments = [self.infile, self.outfile]
        arguments.extend(self.arguments)
        try:
            parsed_arguments = self.parser.parse_args(arguments)
        except SystemExit:
            self.fail("Couldn't parse arguments")
        functions = [f.func for f in parsed_arguments.transforms]
        self.assertEqual(self.functions, functions)

class OrderRespectedTestCase(PopulateTransformsMixIn, unittest.TestCase):
    """
    Ensure that order of arguments translates to order of functions to apply.
    """
    arguments = ['--upper', '--translate', 'dna2protein', '--lower',
            '--squeeze']
    functions = [transform.upper_sequences, transform.translate,
            transform.lower_sequences, transform.squeeze]

class SequenceModTransformsTestCase(PopulateTransformsMixIn, unittest.TestCase):
    arguments = ['--dash-gap',
            '--lower',
            '--reverse',
            '--reverse-complement',
            '--transcribe', 'dna2rna',
            '--translate', 'dna2protein',
            '--ungap',
            '--upper',]
    functions = [transform.dashes_cleanup,
            transform.lower_sequences,
            transform.reverse_sequences,
            transform.reverse_complement_sequences,
            transform.transcribe,
            transform.translate,
            transform.ungap_sequences,
            transform.upper_sequences]

class SeqSelectTransformsTestCase(PopulateTransformsMixIn, unittest.TestCase):

    def setUp(self):
        with tempfile.NamedTemporaryFile(delete=False) as tf:
            self.exclude_from = tf.name
        self.arguments = ['--deduplicate-taxa',
                '--exclude-from-file', self.exclude_from,
                '--include-from-file', self.exclude_from,
                '--head', '10',
                '--max-length', '50',
                '--min-length', '50',
                '--min-ungapped-length', '50',
                '--pattern-include', 'pattern',
                '--pattern-exclude', 'pattern',
                '--prune-empty',
                '--seq-pattern-include', 'pattern',
                '--seq-pattern-exclude', 'pattern',
                ]
        self.functions = [transform.deduplicate_taxa,
                transform.exclude_from_file,
                transform.include_from_file,
                transform.head,
                transform.max_length_discard,
                transform.min_length_discard,
                transform.min_ungap_length_discard,
                transform.name_include,
                transform.name_exclude,
                transform.prune_empty,
                transform.seq_include,
                transform.seq_exclude,
                ]
        super(SeqSelectTransformsTestCase, self).setUp()

    def tearDown(self):
        super(SeqSelectTransformsTestCase, self).tearDown()
        os.remove(self.exclude_from)

class IdModificationTransformsTestCase(PopulateTransformsMixIn, unittest.TestCase):
    arguments = ['--first-name',
            '--name-suffix', 'suffix',
            '--name-prefix', 'prefix',
            '--pattern-replace', '.', 'N',
            '--strip-range']
    functions = [
            transform.first_name_capture,
            transform.name_append_suffix,
            transform.name_insert_prefix,
            transform.name_replace,
            transform.strip_range]

class ArgumentTypeTestCase(PopulateTransformsMixIn, unittest.TestCase):
    arguments = ['--cut', '1:5']
    functions = [transform.multi_cut_sequences]
    def test_argument_type(self):
        arguments = [self.infile, self.outfile]
        arguments.extend(self.arguments)
        try:
            parsed_arguments = self.parser.parse_args(arguments)
        except SystemExit:
            self.fail("Couldn't parse arguments")
        keywords = [f.keywords for f in parsed_arguments.transforms]
        self.assertEqual([{'slices': [slice(0, 5)]}], keywords)

