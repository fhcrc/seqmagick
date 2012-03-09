"""
Filter reads based on quality scores
"""

import argparse
import collections
import csv
import itertools
import re
import sys
from Queue import Queue, Empty
import threading

from Bio import SeqIO
from Bio.SeqIO import QualityIO

from seqmagick import fileformat
from .common import typed_range

# Default minimummean quality score
DEFAULT_MEAN_SCORE = 25.0

# Tools for working with ambiguous bases
# Map from Ambiguous Base to regex
_AMBIGUOUS_MAP = {
       'R': '[GA]',
       'Y': '[TC]',
       'K': '[GT]',
       'M': '[AC]',
       'S': '[GC]',
       'W': '[AT]',
       'B': '[GTC]',
       'D': '[GAT]',
       'H': '[ACT]',
       'V': '[GCA]',
       'N': '[AGCT]',
}

def _ambiguous_pattern(sequence_str):
    """
    Convert an IUPAC ambiguous string to a regular expression pattern string
    """
    return ''.join(_AMBIGUOUS_MAP.get(c, c) for c in sequence_str)


def build_parser(parser):
    """
    Generate a subparser
    """
    parser.add_argument('input_fastq', type=argparse.FileType('r'),
            help="""Input fastq file. A fasta-format file may also be provided
            if --input-qual is also specified.""")
    parser.add_argument('--input-qual', type=argparse.FileType('r'),
            help="""The quality scores associated with the input file. Only
            used if input file is fasta.""")
    parser.add_argument('output_file', type=argparse.FileType('w'),
            help="""Output file. Format determined from extension.""")
    output_group = parser.add_argument_group("Output")

    output_group.add_argument('--report-out', type=argparse.FileType('w'),
            default=sys.stdout, help="""Output file for report [default:
            stdout]""")
    output_group.add_argument('--failure-out', type=argparse.FileType('w'),
            help="""File to write failure report [default: None]""")

    parser.add_argument('--min-mean-quality', metavar='QUALITY', type=float,
            default=DEFAULT_MEAN_SCORE, help="""Minimum mean quality score for
            each read [default: %(default)s]""")
    parser.add_argument('--min-length', metavar='LENGTH', type=int,
            default=200, help="""Minimum length to keep sequence [default:
            %(default)s]""")
    parser.add_argument('--max-length', metavar='LENGTH', type=int,
            default=1000, help="""Maximum length to keep before truncating
            [default: %(default)s]. This operation occurs before
            --max-ambiguous""")


    window_group = parser.add_argument_group('Quality window options')
    window_group.add_argument('--quality-window-mean-qual', type=float,
            help="""Minimum quality score within the window defined by
            --quality-window. [default: same as --min-mean-quality]""")
    window_group.add_argument('--quality-window-prop', help="""Proportion of
            reads within quality window to that must pass filter. Floats are [default:
            %(default).1f]""", default=1.0, type=typed_range(float, 0.0, 1.0))
    window_group.add_argument('--quality-window', type=int, metavar='WINDOW_SIZE',
            default=0, help="""Window size for truncating sequences.  When set
            to a non-zero value, sequences are truncated where the mean mean
            quality within the window drops below --min-mean-quality.
            [default: %(default)s]""")

    parser.add_argument('--ambiguous-action', choices=('truncate', 'drop'),
            help="""Action to take on ambiguous base in sequence (N's).
            [default: no action]""")
    parser.add_argument('--max-ambiguous', default=None, help="""Maximum number
            of ambiguous bases in a sequence. Sequences exceeding this count
            will be removed.""", type=int)

    barcode_group = parser.add_argument_group('Barcode/Primer')
    barcode_group.add_argument('--primer', help="""IUPAC ambiguous primer to
            require""")
    barcode_group.add_argument('--barcode-file', help="""CSV file
            containing sample_id,barcode rows""", type=argparse.FileType('r'))
    barcode_group.add_argument('--barcode-header', action='store_true',
            default=False, help="""Barcodes have a header row [default:
            %(default)s]""")
    barcode_group.add_argument('--map-out', help="""Path to write
            sequence_id,sample_id pairs""", type=argparse.FileType('w'),
            metavar='SAMPLE_MAP')
    barcode_group.add_argument('--quoting', help="""A string naming an
            attribute of the csv module defining the quoting behavior for
            `SAMPLE_MAP`.  [default: %(default)s]""", default='QUOTE_MINIMAL',
            choices=[s for s in dir(csv) if s.startswith('QUOTE_')])

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

class FailureReportWriter(threading.Thread):
    """
    Writes a log of sequences that failed filtering, and the filter that
    removed them.
    """
    def __init__(self, queue, fp):
        super(FailureReportWriter, self).__init__()
        self.queue = queue
        self.writer = csv.DictWriter(fp, ('failed_sequence', 'reason', 'value'),
                delimiter='\t', lineterminator='\n')
        self.writer.writeheader()

    def run(self):
        while True:
            try:
                record = self.queue.get(timeout=1)
            except Empty:
                continue

            if record is None:
                self.queue.task_done()
                return
            self.writer.writerow(record)
            self.queue.task_done()


class Failure(object):
    """
    Simple object to hold a value corresponding to a failure.  Evaluates to
    false in boolean tests.
    """
    def __init__(self, value=None):
        self.value = value

    def __nonzero__(self):
        return False

class BaseFilter(object):
    """
    Base class for filters
    """
    report_fields = ['name', 'passed_unchanged', 'passed_changed', 'failed',
            'total_filtered', 'proportion_passed']

    def __init__(self):
        self.passed_unchanged = 0
        self.passed_changed = 0
        self.failed = 0

    def filter_record(self, record):
        """
        Filter a record. If the filter succeeds, returns a SeqRecord. If it
        fails, returns either None or an instance of Failure containing the
        value with failed.
        """
        raise NotImplementedError("Override in subclass")

    def filter_records(self, records, failure_queue=None):
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
                if failure_queue:
                    value = filtered.value if filtered is not None else None
                    failure_queue.put({'failed_sequence': record.id,
                        'reason': self.name, 'value': value})

    @property
    def passed(self):
        return self.passed_changed + self.passed_unchanged

    @property
    def total_filtered(self):
        return self.passed + self.failed

    @property
    def proportion_passed(self):
        if not self.total_filtered:
            return 0
        return float(self.passed) / self.total_filtered

    def report_dict(self):
        return dict((f, getattr(self, f)) for f in self.report_fields)


class QualityScoreFilter(BaseFilter):
    """
    Quality score filter - requires that the average base quality over the
    length of the read is greater than some threshold.
    """

    def __init__(self, min_mean_score=DEFAULT_MEAN_SCORE):
        super(QualityScoreFilter, self).__init__()
        self.min_mean_score = min_mean_score
        self.name = "Quality Score [min_mean: {0}]".format(min_mean_score)

    def filter_record(self, record):
        """
        Filter a single record

        Returns None if the record failed.
        """
        quality_scores = record.letter_annotations['phred_quality']

        mean_score = mean(quality_scores)
        return record if mean_score >= self.min_mean_score else Failure(mean_score)

class WindowQualityScoreFilter(BaseFilter):
    """
    Filter records, truncating records when the mean score drops below a
    certain value.
    """
    def __init__(self, window_size, min_mean_score=DEFAULT_MEAN_SCORE):
        super(WindowQualityScoreFilter, self).__init__()
        self.min_mean_score = min_mean_score
        assert window_size and window_size > 0
        self.window_size = window_size
        self.name = ("Windowed Quality Score " +
                     "[min_mean-quality: {0}; window_size: {1}]").format(
                             min_mean_score, window_size)

    def filter_record(self, record):
        """
        Filter a single record

        Returns None if the record failed.
        """
        quality_scores = record.letter_annotations['phred_quality']

        # Simple case - window covers whole sequence
        if len(record) <= self.window_size:
            mean_score = mean(quality_scores)
            return record if mean_score >= self.min_mean_score else Failure(mean_score)

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

        if clip_right:
            return record[:clip_right]

class AmbiguousBaseFilter(BaseFilter):
    """
    Filter records, taking some action if 'N' is encountered in the sequence.

    action  - either 'truncate' (drop N and any sequence following) or 'drop'
              (remove sequences with 'N's)
    """
    name = 'Ambiguous Base'

    def __init__(self, action):
        super(AmbiguousBaseFilter, self).__init__()
        if action not in ('truncate', 'drop'):
            raise ValueError("Unknown action: {0}".format(action))
        self.action = action
        self.name = AmbiguousBaseFilter.name + " [{0}]".format(action)

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

class MaxAmbiguousFilter(BaseFilter):
    """
    Filters records exceeding some minimum number of ambiguous bases
    """
    name = "Maximum Ambiguous Bases"

    def __init__(self, max_ambiguous):
        super(MaxAmbiguousFilter, self).__init__()
        assert max_ambiguous is not None
        self.max_ambiguous = max_ambiguous
        self.name = self.name + '[{0}]'.format(max_ambiguous)

    def filter_record(self, record):
        n_count = record.seq.upper().count('N')
        if n_count > self.max_ambiguous:
            return Failure(n_count)
        else:
            assert n_count <= self.max_ambiguous
            return record

class MinLengthFilter(BaseFilter):
    """
    Remove records which don't meet minimum length
    """
    def __init__(self, min_length):
        super(MinLengthFilter, self).__init__()
        assert min_length > 0
        self.min_length = min_length
        self.name = "Minimum Length [{0}]".format(min_length)

    def filter_record(self, record):
        """
        Filter record, dropping any that don't meet minimum length
        """
        l = len(record)
        if l >= self.min_length:
            return record
        else:
            return Failure(l)

class MaxLengthFilter(BaseFilter):
    """
    Truncate long sequences
    """
    name = "Maximum Length"
    def __init__(self, max_length):
        super(MaxLengthFilter, self).__init__()
        self.max_length = max_length
        self.name = self.name + " [{0}]".format(max_length)

    def filter_record(self, record):
        """
        Filter record, truncating any over some maximum length
        """
        if len(record) >= self.max_length:
            return record[:self.max_length]
        else:
            return record

class PrimerBarcodeFilter(BaseFilter):
    """
    Filter that checks that the sequence starts with a known barcode/primer
    combination.

    Sequences that pass the filter have the barcode and primer removed.

    If an output_file is provided, (sequence_id, sample_id) tuples are written
    to it.
    """
    name = "Primer/Barcode"

    def __init__(self, primer, barcodes, output_file=None, trim=True, quoting=csv.QUOTE_MINIMAL):
        super(PrimerBarcodeFilter, self).__init__()
        self.primer = primer
        self.barcodes = barcodes
        self.trim = True
        if output_file:
            self.writer = csv.writer(output_file, lineterminator='\n', quoting = quoting)
        else:
            self.writer = None

        self.pattern = re.compile('^({0}){1}'.format('|'.join(barcodes),
                                                    _ambiguous_pattern(primer)),
                                  re.IGNORECASE)

    def _report_match(self, record, sample):
        if not self.writer:
            return
        self.writer.writerow((record.id, sample))

    def filter_record(self, record):
        m = self.pattern.match(str(record.seq))
        if m:
            self._report_match(record, self.barcodes[m.group(1)])
            if self.trim:
                record = record[m.end():]
            return record

def parse_barcode_file(fp, header=False):
    """
    Load label, barcode pairs from a CSV file.

    Returns a map from barcode -> label

    Any additional columns are ignored
    """
    d = {}
    reader = csv.reader(fp)

    if header:
        # Skip header
        next(reader)

    # Skip blank rows
    records = (record for record in reader if record)

    for record in records:
        value, key = record[:2]
        if key in d:
            raise ValueError("Duplicate barcode: {0}".format(key))
        if value in d.values():
            raise ValueError("Duplicate sample label: {0}".format(value))
        d[key] = value

    # Require that all barcodes are the same length
    if min(len(k) for k in d) != max(len(k) for k in d):
        keys = d.keys()
        keys.sort(key=len)
        minlen_key = keys[0]
        maxlen_key = keys[-1]

        raise ValueError(("Length of barcode '{0}' does not match "
            "length of barcode '{1}'").format(minlen_key, maxlen_key))

    return d

def action(arguments):
    """
    Given parsed arguments, filter input files.
    """
    if arguments.quality_window_mean_qual and not arguments.quality_window:
        raise ValueError("--quality-window-mean-qual specified without "
                "--quality-window")

    queue = None
    if arguments.failure_out:
        queue = Queue()
        t = FailureReportWriter(queue, arguments.failure_out)
        t.setDaemon(True)
        t.start()

    # Always filter with a quality score
    qfilter = QualityScoreFilter(arguments.min_mean_quality)
    filters = [qfilter]

    output_type = fileformat.from_filename(arguments.output_file.name)
    with arguments.input_fastq as fp:
        if arguments.input_qual:
            sequences = QualityIO.PairedFastaQualIterator(fp,
                    arguments.input_qual)
        else:
            sequences = SeqIO.parse(fp, 'fastq')

        # Add filters
        if arguments.max_length:
            max_length_filter = MaxLengthFilter(arguments.max_length)
            filters.append(max_length_filter)
        if arguments.min_length:
            min_length_filter = MinLengthFilter(arguments.min_length)
            filters.append(min_length_filter)
        if arguments.max_ambiguous is not None:
            max_ambig_filter = MaxAmbiguousFilter(arguments.max_ambiguous)
            filters.append(max_ambig_filter)
        if arguments.ambiguous_action:
            ambiguous_filter = AmbiguousBaseFilter(
                    arguments.ambiguous_action)
            filters.append(ambiguous_filter)
        if arguments.quality_window:
            min_qual = arguments.quality_window_mean_qual or \
                    arguments.min_mean_quality
            window_filter = WindowQualityScoreFilter(arguments.quality_window,
                    min_qual)
            filters.insert(0, window_filter)

        if arguments.barcode_file:
            with arguments.barcode_file:
                barcodes = parse_barcode_file(arguments.barcode_file,
                        arguments.barcode_header)
            f = PrimerBarcodeFilter(arguments.primer or '', barcodes,
                    arguments.map_out,
                    quoting=getattr(csv, arguments.quoting))
            filters.append(f)

        for f in filters:
            sequences = f.filter_records(sequences, queue)

        with arguments.output_file:
            SeqIO.write(sequences, arguments.output_file, output_type)

    rpt_rows = (f.report_dict() for f in filters)

    # Write report
    with arguments.report_out as fp:
        writer = csv.DictWriter(fp, BaseFilter.report_fields,
                lineterminator='\n', delimiter='\t')
        writer.writeheader()
        writer.writerows(rpt_rows)

    if queue:
        queue.join()
