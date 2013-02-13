"""
Filter reads based on quality scores
"""

import collections
import csv
import itertools
import logging
import os
import sys
import time

from Bio import SeqIO
try:
    from Bio import trie, triefind
except ImportError:
    trie = None
    triefind = None
from Bio.SeqIO import QualityIO

from seqmagick import fileformat, __version__
from .common import typed_range, FileType

# Default minimummean quality score
DEFAULT_MEAN_SCORE = 25.0

# Tools for working with ambiguous bases
# Map from Ambiguous Base to regex
_AMBIGUOUS_MAP = {
       'R': 'GA',
       'Y': 'TC',
       'K': 'GT',
       'M': 'AC',
       'S': 'GC',
       'W': 'AT',
       'B': 'GTC',
       'D': 'GAT',
       'H': 'ACT',
       'V': 'GCA',
       'N': 'AGCT',
}

def all_unambiguous(sequence_str):
    """
    All unambiguous versions of sequence_str
    """
    result = [[]]
    for c in sequence_str:
        result = [i + [a] for i in result
                  for a in _AMBIGUOUS_MAP.get(c, c)]
    return [''.join(i) for i in result]

def build_parser(parser):
    """
    Generate a subparser
    """
    parser.add_argument('input_fastq', type=FileType('r'),
            help="""Input fastq file. A fasta-format file may also be provided
            if --input-qual is also specified.""")
    parser.add_argument('--input-qual', type=FileType('r'),
            help="""The quality scores associated with the input file. Only
            used if input file is fasta.""")
    parser.add_argument('output_file', type=FileType('w'),
            help="""Output file. Format determined from extension.""")

    output_group = parser.add_argument_group("Output")
    output_group.add_argument('--report-out', type=FileType('w'),
            default=sys.stdout, help="""Output file for report [default:
            stdout]""")
    output_group.add_argument('--details-out', type=FileType('w'),
             help="""Output file to report fate of each sequence""")
    output_group.add_argument('--no-details-comment', action='store_false',
            default=True, dest='details_comment', help="""Do not write comment
            lines with version and call to start --details-out""")

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
    primer_group = barcode_group.add_mutually_exclusive_group()
    primer_group.add_argument('--primer', help="""IUPAC ambiguous primer to
            require""")
    primer_group.add_argument('--no-primer', help="""Do not use a primer.""",
            action='store_const', const='', dest='primer')
    barcode_group.add_argument('--barcode-file', help="""CSV file containing
            sample_id,barcode[,primer] in the rows. A single primer for all
            sequences may be specified with `--primer`, or `--no-primer` may be
            used to indicate barcodes should be used without a primer
            check.""", type=FileType('r'))
    barcode_group.add_argument('--barcode-header', action='store_true',
            default=False, help="""Barcodes have a header row [default:
            %(default)s]""")
    barcode_group.add_argument('--map-out', help="""Path to write
            sequence_id,sample_id pairs""", type=FileType('w'),
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

class FailedFilter(Exception):
    """
    A read failed filtering
    """
    def __init__(self, value=None):
        self.value = value

class RecordEventListener(object):
    """
    Contains and dispatches to handlers on events around sequence records

    Event handlers take a single positional argument, the record, and optional
    additional keyword arguments.
    """
    def __init__(self):
        self.listeners = collections.defaultdict(set)

    def __call__(self, event, record, **kwargs):
        """
        Trigger an event

        :param event: Event name
        :param record: Record affected
        :param **kwargs: Optional additional arguments to pass to handlers
        """
        if event in self.listeners:
            for listener in self.listeners[event]:
                listener(record, **kwargs)

    def register_handler(self, event, handler):
        """
        Register ``handler`` for ``event``
        """
        self.listeners[event].add(handler)

    def iterable_hook(self, name, iterable):
        """
        Fire an event named ``name`` with each item in iterable
        """
        for record in iterable:
            self(name, record)
            yield record

class RecordReportHandler(object):
    """
    Generates a report to a CSV file detailing every record processed.

    Listens for events: [read, write, failed_filter, found_barcode]
    """
    HEADERS = ('sequence_name', 'in_length', 'in_mean_qual', 'sample',
               'out_length', 'out_mean_qual', 'fail_filter', 'fail_value')
    def __init__(self, fp, args, write_comments=True):
        if write_comments:
            fp.write('# Generated by `seqmagick quality-filter` version {0}\n'.format(__version__))
            fp.write('# Arguments: {0}\n'.format(' '.join(args)))
            fp.write('# Working directory: {0}\n'.format(os.getcwd()))

        self.writer = csv.DictWriter(fp, self.HEADERS, lineterminator='\n',
                quoting=csv.QUOTE_NONNUMERIC)
        self.writer.writeheader()
        self.current_record = None

        self.read = 0
        self.failed = 0
        self.start = time.time()
        self.last_report = 0.0

    def register_with(self, listener):
        listener.register_handler('failed_filter', self._record_failed)
        listener.register_handler('read', self._read_record)
        listener.register_handler('write', self._wrote_record)
        listener.register_handler('found_barcode', self._found_barcode)

    def _write(self):
        assert self.current_record
        self.writer.writerow(self.current_record)
        self.current_record = None

    def _record_failed(self, record, filter_name, value=None):
        self.current_record.update({'fail_filter': filter_name,
                                    'fail_value': value})

        self._write()
        self.failed += 1
        self._report()

    def _read_record(self, record):
        self.current_record = {'sequence_name': record.id,
                               'in_length': len(record)}
        if 'phred_quality' in record.letter_annotations:
            self.current_record['in_mean_qual'] = \
                    mean(record.letter_annotations['phred_quality'])
        self.read += 1

    def _found_barcode(self, record, sample, barcode=None):
        """Hook called when barcode is found"""
        assert record.id == self.current_record['sequence_name']
        self.current_record['sample'] = sample

    def _wrote_record(self, record):
        self.current_record['out_length'] = len(record)
        if 'phred_quality' in record.letter_annotations:
            self.current_record['out_mean_qual'] = \
                    mean(record.letter_annotations['phred_quality'])
        self._write()
        self._report()

    def _report(self):
        if not sys.stdout.isatty():
            return
        t = time.time()
        if t - self.last_report < 0.4 or not self.read:
            return

        self.last_report = t
        sys.stderr.write('{0:10.1f}s Processed {1:10d} records; {2:10d} passed ({3:6.2f}%)\r'.format(
            t - self.start,
            self.read,
            self.read - self.failed,
            float(self.read - self.failed) / self.read * 100.0))


class BaseFilter(object):
    """
    Base class for filters
    """
    report_fields = ('name', 'passed_unchanged', 'passed_changed', 'failed',
                     'total_filtered', 'proportion_passed')

    def __init__(self, listener=None):
        self.passed_unchanged = 0
        self.passed_changed = 0
        self.failed = 0
        self.listener = listener

    def filter_record(self, record):
        """
        Filter a record. If the filter succeeds, returns a SeqRecord. If it
        fails, raises an instance of FailedFilter with an optional value.
        """
        raise NotImplementedError("Override in subclass")

    def filter_records(self, records):
        """
        Apply the filter to records
        """
        for record in records:
            try:
                filtered = self.filter_record(record)
                assert(filtered)
                # Quick tracking whether the sequence was modified
                if filtered == record:
                    self.passed_unchanged += 1
                else:
                    self.passed_changed += 1
                yield filtered
            except FailedFilter as e:
                self.failed += 1
                v = e.value
                if self.listener:
                    self.listener('failed_filter', record, filter_name=self.name, value=v)

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
        """
        quality_scores = record.letter_annotations['phred_quality']

        mean_score = mean(quality_scores)
        if mean_score >= self.min_mean_score:
            return record
        else:
            raise FailedFilter(mean_score)

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
        """
        quality_scores = record.letter_annotations['phred_quality']

        # Simple case - window covers whole sequence
        if len(record) <= self.window_size:
            mean_score = mean(quality_scores)
            if mean_score >= self.min_mean_score:
                return record
            else:
                raise FailedFilter(mean_score)

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
        else:
            # First window failed - record fails
            raise FailedFilter()

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
            raise FailedFilter()
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
            raise FailedFilter(n_count)
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
            raise FailedFilter(l)

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

    def __init__(self, trie, output_file=None, trim=True, quoting=csv.QUOTE_MINIMAL):
        super(PrimerBarcodeFilter, self).__init__()
        self.trim = True
        self.trie = trie

    def filter_record(self, record):
        m = triefind.match(str(record.seq), self.trie)
        if m:
            if self.listener:
                self.listener('found_barcode', record, barcode=m, sample=self.trie[m])
            if self.trim:
                record = record[len(m):]
            return record
        else:
            raise FailedFilter()

def parse_barcode_file(fp, primer=None, header=False):
    """
    Load label, barcode, primer records from a CSV file.

    Returns a map from barcode -> label

    Any additional columns are ignored
    """
    tr = trie.trie()
    reader = csv.reader(fp)

    if header:
        # Skip header
        next(reader)

    # Skip blank rows
    records = (record for record in reader if record)

    for record in records:
        specimen, barcode = record[:2]
        if primer is not None:
            pr = primer
        else:
            pr = record[2]
        for sequence in all_unambiguous(barcode + pr):
            if tr.has_key(sequence):
                raise ValueError("Duplicate sample: {0}, {1} both have {2}",
                        specimen, tr[sequence], sequence)
            logging.info('%s->%s', sequence, specimen)
            tr[sequence] = specimen

    return tr

def action(arguments):
    """
    Given parsed arguments, filter input files.
    """
    if arguments.quality_window_mean_qual and not arguments.quality_window:
        raise ValueError("--quality-window-mean-qual specified without "
                "--quality-window")

    if trie is None or triefind is None:
        raise ValueError('Missing Bio.trie and/or Bio.triefind modules. Cannot continue')

    # Always filter with a quality score
    qfilter = QualityScoreFilter(arguments.min_mean_quality)
    filters = [qfilter]

    output_type = fileformat.from_handle(arguments.output_file)
    with arguments.input_fastq as fp:
        if arguments.input_qual:
            sequences = QualityIO.PairedFastaQualIterator(fp,
                    arguments.input_qual)
        else:
            sequences = SeqIO.parse(fp, 'fastq')

        listener = RecordEventListener()
        if arguments.details_out:
            rh = RecordReportHandler(arguments.details_out, arguments.argv,
                    arguments.details_comment)
            rh.register_with(listener)

        # Track read sequences
        sequences = listener.iterable_hook('read', sequences)

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
                tr = parse_barcode_file(arguments.barcode_file,
                        arguments.primer, arguments.barcode_header)
            f = PrimerBarcodeFilter(tr)
            filters.append(f)

            if arguments.map_out:
                barcode_writer = csv.writer(arguments.map_out,
                        quoting=getattr(csv, arguments.quoting),
                        lineterminator='\n')
                def barcode_handler(record, sample, barcode=None):
                    barcode_writer.writerow((record.id, sample))
                listener.register_handler('found_barcode', barcode_handler)
        for f in filters:
            f.listener = listener
            sequences = f.filter_records(sequences)

        # Track sequences which passed all filters
        sequences = listener.iterable_hook('write', sequences)

        with arguments.output_file:
            SeqIO.write(sequences, arguments.output_file, output_type)

    rpt_rows = (f.report_dict() for f in filters)

    # Write report
    with arguments.report_out as fp:
        writer = csv.DictWriter(fp, BaseFilter.report_fields,
                lineterminator='\n', delimiter='\t')
        writer.writeheader()
        writer.writerows(rpt_rows)
