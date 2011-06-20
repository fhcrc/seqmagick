"""
Filter reads based on quality scores
"""

import argparse
import collections
import csv
import itertools
import sys

from Bio import SeqIO
from Bio.SeqIO import QualityIO

from seqmagick import fileformat

def build_parser(parser):
    """
    Generate a subparser
    """
    parser.add_argument('input_fasta', type=argparse.FileType('r'),
            help='Input fasta file')
    parser.add_argument('input_qual', type=argparse.FileType('r'),
            help='The quality scores associated with fasta_file')
    parser.add_argument('output_file', type=argparse.FileType('w'),
            help="""Output file. Format determined from extension.""")
    parser.add_argument('--min-mean-quality', metavar='QUALITY', type=float,
            default=25, help="""Minimum mean quality score for each read
            [default: %(default)s]""")
    parser.add_argument('--min-length', metavar='LENGTH', type=int,
            help="""Minimum length to keep sequence [default: %(default)s]""")
    parser.add_argument('--quality-window', type=int, metavar='WINDOW_SIZE',
            default=0, help="""Window size for truncating sequences.  When set
            to a non-zero value, sequences are truncated where the mean mean
            quality within the window drops below --min-mean-quality.
            [default: %(default)s]""")

    parser.add_argument('--ambiguous-action', choices=('truncate', 'drop'),
            help="""Action to take on ambiguous base in sequence (N's).
            [default: no action]""")


def mean(sequence):
    """
    Calculates the arithmetic mean of a list / tuple
    """
    return sum(sequence) / float(len(sequence))


def moving_average(iterable, n):
    """
    From Python collections module documentation

    moving_average([40, 30, 50, 46, 39, 44]) --> 40.0 42.0 45.0 43.0
    """
    it = iter(iterable)
    d = collections.deque(itertools.islice(it, n-1))
    d.appendleft(0)
    s = sum(d)
    for elem in it:
        s += elem - d.popleft()
        d.append(elem)
        yield s / float(n)


class BaseFilter(object):
    """
    Base class for filters
    """

    def __init__(self):
        self.passed_unchanged = 0
        self.passed_changed = 0
        self.failed = 0

    def filter_record(self, record):
        raise NotImplementedError("Override in subclass")

    def filter_records(self, records):
        """
        Apply the filter to records
        """
        for record in records:
            filtered = self.filter_record(record)
            if filtered:
                # Quick tracking whether the sequence was modified
                if filtered == record:
                    self.passed_unchanged += 1
                else:
                    self.passed_changed += 1
                yield filtered
            else:
                self.failed += 1

    @property
    def passed(self):
        return self.passed_changed + self.passed_unchanged

    @property
    def total_filtered(self):
        return self.passed + self.failed

    @property
    def proportion_passed(self):
        return float(self.passed) / self.total_filtered


class QualityScoreFilter(BaseFilter):
    """
    Quality score filter
    """

    name = 'Quality Score Filter'

    def __init__(self, min_mean_score=25.0, window_size=0):
        super(QualityScoreFilter, self).__init__()
        self.min_mean_score = min_mean_score
        self.window_size = window_size

    def filter_record(self, record):
        """
        Filter a single record

        Returns None if the record failed.
        """
        quality_scores = record.letter_annotations['phred_quality']

        # Simple case - no window size or window covers whole sequence
        if self.window_size == 0 or len(record) <= self.window_size:
            mean_score = mean(quality_scores)
            return record if mean_score >= self.min_mean_score else None

        # Find the right clipping point. Start clipping at the beginning of the
        # sequence, then extend the window to include regions with acceptable
        # mean quality scores.
        clip_right = 0
        for i, a in enumerate(moving_average(quality_scores,
            self.window_size)):

            if a >= self.min_mean_score:
                clip_right = i + self.window_size
            else:
                break

        return record[:clip_right]


class AmbiguousBaseFilter(BaseFilter):
    """
    Filter records, taking some action if 'N' is encountered in the sequence.

    action  - either 'truncate' (drop N and any sequence following) or 'drop'
              (remove sequences with 'N's)
    """

    name = 'Ambiguous Base Filter'

    def __init__(self, action):
        super(AmbiguousBaseFilter, self).__init__()
        if action not in ('truncate', 'drop'):
            raise ValueError("Unknown action: {0}".format(action))
        self.action = action

    def filter_record(self, record):
        """
        Filter a record, truncating or dropping at an 'N'
        """
        nloc = record.seq.find('N')
        if nloc == -1:
            return record
        elif self.action == 'truncate':
            return record[:nloc]
        elif self.action == 'drop':
            return None
        else:
            assert False

class MinLengthFilter(BaseFilter):
    def __init__(self, min_length):
        self.min_length = min_length

    def filter_record(self, record):
        """
        Filter record, dropping any that don't meet minimum lenght
        """
        if len(record) >= self.min_length:
            return record

def action(arguments):
    """
    Given parsed arguments, filter input files.
    """
    qfilter = QualityScoreFilter(arguments.min_mean_quality,
            arguments.quality_window)

    output_type = fileformat.from_filename(arguments.output_file.name)
    filters = [qfilter]
    with arguments.input_fasta:
        with arguments.input_qual:
            sequences = QualityIO.PairedFastaQualIterator(
                    arguments.input_fasta, arguments.input_qual)
            filtered = qfilter.filter_records(sequences)
            if arguments.ambiguous_action:
                ambiguous_filter = AmbiguousBaseFilter(
                        arguments.ambiguous_action)
                filtered = ambiguous_filter.filter_records(filtered)
                filters.append(ambiguous_filter)
            if arguments.min_length:
                min_length_filter = MinLengthFilter(arguments.min_length)
                filtered = min_length_filter.filter_records(filtered)
                filters.append(min_length_filter)

            with arguments.output_file:
                SeqIO.write(filtered, arguments.output_file,
                        output_type)

    rpt_rows = [(f.name, f.passed_unchanged, f.passed_changed, f.failed,
        f.total_filtered, f.proportion_passed) for f in filters]

    # Write report
    writer = csv.writer(sys.stdout, lineterminator='\n', delimiter='\t')
    writer.writerow(('filter', 'passed_unchanged', 'passed_changed', 'failed',
        'total_processed', 'proportion_passed'))
    writer.writerows(rpt_rows)
