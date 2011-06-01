"""
Filter reads based on quality scores
"""

import argparse
import collections
import itertools

from Bio import SeqIO
from Bio.SeqIO import QualityIO


def build_parser(parser):
    """
    Generate a subparser
    """
    parser.add_argument('input_fasta', type=argparse.FileType('r'),
            help='Input fasta file')
    parser.add_argument('input_qual', type=argparse.FileType('r'),
            help='The quality scores associated with fasta_file')
    parser.add_argument('output_fasta', type=argparse.FileType('w'),
            help='Output FASTA')

    parser.add_argument('--min-mean-quality', type=float,
            default=25.0, help="""Minimum mean quality score for each read
            (default: %(default)s""")
    parser.add_argument('--quality-window', type=int,
            default=0, help="""Window size for truncating sequences once the
            mean quality within the window drops below --min-mean-quality.
            0=use whole sequence. (default: %(default)s)""")

    parser.add_argument('--ambiguous-action', choices=('truncate', 'drop'),
            help="""Action to take on ambiguous base in sequence (N's).
            Default: no action.""")

def mean(sequence):
    return sum(sequence) / float(len(sequence))

def moving_average(iterable, n=3):
    # moving_average([40, 30, 50, 46, 39, 44]) --> 40.0 42.0 45.0 43.0
    # http://en.wikipedia.org/wiki/Moving_average
    it = iter(iterable)
    d = collections.deque(itertools.islice(it, n-1))
    d.appendleft(0)
    s = sum(d)
    for elem in it:
        s += elem - d.popleft()
        d.append(elem)
        yield s / float(n)


class QualityScoreFilter(object):
    """
    Quality score filter
    """

    def __init__(self, min_mean_score=25.0, window_size=0):
        self.min_mean_score = min_mean_score
        self.window_size = window_size

    def filter_record(self, record):
        """
        Filter a single record

        Returns None
        """
        quality_scores = record.letter_annotations['phred_quality']
        # Simple case - no window size or window covers whole sequence
        if self.window_size == 0 or len(record) <= self.window_size:
            mean_score = mean(quality_scores)
            return record if mean_score >= self.min_mean_score else None

        # Find the right clipping point
        clip_right = 0
        for i, a in enumerate(moving_average(quality_scores,
            self.window_size)):

            if a >= self.min_mean_score:
                clip_right = i + self.window_size
            else:
                break

        return record[:clip_right]

    def filter_records(self, records):
        """
        Apply the filter to records
        """
        for record in records:
            filtered = self.filter_record(record)
            if filtered:
                yield filtered


def ambiguous_base_filter(records, action):
    """
    Filter records, taking some action if 'N' is encountered in the sequence.

    records - iterable of SeqRecord objects
    action  - either 'truncate' (drop N and any sequence following) or 'drop'
              (remove sequences with 'N's)
    """
    if action not in ('truncate', 'drop'):
        raise ValueError("Unknown action: {0}".format(action))
    for record in records:
        nloc = record.seq.find('N')
        if nloc == -1:
            yield record
        elif action == 'truncate':
            yield record[:nloc]
        elif action == 'drop':
            continue
        else:
            assert False

def action(arguments):
    qfilter = QualityScoreFilter(arguments.min_mean_quality,
            arguments.quality_window)
    with arguments.input_fasta:
        with arguments.input_qual:
            sequences = QualityIO.PairedFastaQualIterator(
                    arguments.input_fasta, arguments.input_qual)
            filtered = qfilter.filter_records(sequences)
            if arguments.ambiguous_action:
                filtered = ambiguous_base_filter(filtered,
                        arguments.ambiguous_action)

            with arguments.output_fasta:
                count = SeqIO.write(filtered, arguments.output_fasta,
                        'fasta')

            print count
