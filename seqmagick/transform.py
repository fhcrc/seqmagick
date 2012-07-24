"""
Functions to transform / filter sequences
"""
import collections
import contextlib
import cPickle as pickle
import itertools
import logging
import re
import string
import tempfile

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid

# Characters to be treated as gaps
GAP_CHARS = "-."

# Size of temporary file buffer: default to 20MB
DEFAULT_BUFFER_SIZE = 20971520 # 20*2**20

@contextlib.contextmanager
def _record_buffer(records, buffer_size=DEFAULT_BUFFER_SIZE):
    """
    Buffer for transform functions which require multiple passes through data.

    Value returned by context manager is a function which returns an iterator
    through records.
    """
    with tempfile.SpooledTemporaryFile(buffer_size, mode='wb+') as tf:
        pickler = pickle.Pickler(tf)
        for record in records:
            pickler.dump(record)

        def record_iter():
            tf.seek(0)
            unpickler = pickle.Unpickler(tf)
            while True:
                try:
                    yield unpickler.load()
                except EOFError:
                    break

        yield record_iter


def dashes_cleanup(records, prune_chars='.:?~'):
    """
    Take an alignment and convert any undesirable characters such as ? or ~ to
    -.
    """
    logging.info("Applying _dashes_cleanup: converting any . or : to -.")
    translation_table = string.maketrans(prune_chars, '-' * len(prune_chars))
    for record in records:
        record.seq = Seq(str(record.seq).translate(translation_table),
                         record.seq.alphabet)
        yield record


def deduplicate_sequences(records, out_file):
    """
    Remove any duplicate records with identical sequences, keep the first
    instance seen and discard additional occurences.
    """

    logging.info('Applying _deduplicate_sequences generator: '
                 'removing any duplicate records with identical sequences.')
    checksum_sequences = collections.defaultdict(list)
    for record in records:
        checksum = seguid(record.seq)
        sequences = checksum_sequences[checksum]
        if not sequences:
            yield record
        sequences.append(record.id)

    if out_file is not None:
        with out_file:
            for sequences in checksum_sequences.itervalues():
                out_file.write('%s\n' % (' '.join(sequences),))


def deduplicate_taxa(records):
    """
    Remove any duplicate records with identical IDs, keep the first
    instance seen and discard additional occurences.
    """
    logging.info('Applying _deduplicate_taxa generator: ' + \
                 'removing any duplicate records with identical IDs.')
    taxa = set()
    for record in records:
        # Default to full ID, split if | is found.
        taxid = record.id
        if '|' in record.id:
            try:
                taxid = int(record.id.split("|")[0])
            except:
                # If we couldn't parse an integer from the ID, just fall back
                # on the ID
                logging.warn("Unable to parse integer taxid from %s",
                        taxid)
        if taxid in taxa:
            continue
        taxa.add(taxid)
        yield record


def first_name_capture(records):
    """
    Take only the first whitespace-delimited word as the name of the sequence.
    Essentially removes any extra text from the sequence's description.
    """
    logging.info('Applying _first_name_capture generator: '
                 'making sure ID only contains the first whitespace-delimited '
                 'word.')
    whitespace = re.compile(r'\s+')
    for record in records:
        if whitespace.search(record.description):
            yield SeqRecord(record.seq, id=record.id,
                            description="")
        else:
            yield record


def include_from_file(records, handle):
    """
    Filter the records, keeping only sequences whose ID is contained in the
    handle.
    """
    ids = set(i.strip() for i in handle)

    for record in records:
        if record.id.strip() in ids:
            yield record


def exclude_from_file(records, handle):
    """
    Filter the records, keeping only sequences whose ID is not contained in the
    handle.
    """
    ids = set(i.strip() for i in handle)

    for record in records:
        if record.id.strip() not in ids:
            yield record


def isolate_region(sequences, start, end, gap_char='-'):
    """
    Replace regions before and after start:end with gap chars
    """
    # Check arguments
    if end <= start:
        raise ValueError("start of slice must precede end ({0} !> {1})".format(
            end, start))

    for sequence in sequences:
        seq = sequence.seq
        start_gap = gap_char * start
        end_gap = gap_char * (len(seq) - end)
        seq = Seq(start_gap + str(seq[start:end]) + end_gap,
                alphabet=seq.alphabet)
        sequence.seq = seq
        yield sequence


def _cut_sequences(records, cut_slice):
    """
    Cut sequences given a slice.
    """
    for record in records:
        yield record[cut_slice]

def drop_columns(records, slices):
    """
    Drop all columns present in ``slices`` from records
    """
    for record in records:
        # Generate a set of indices to remove
        drop = set(i for slice in slices
                   for i in range(*slice.indices(len(record))))
        keep = [i not in drop for i in xrange(len(record))]
        record.seq = Seq(''.join(itertools.compress(record.seq, keep)), record.seq.alphabet)
        yield record

def multi_cut_sequences(records, slices):
    # If only a single slice is specified, use _cut_sequences,
    # since this preserves per-letter annotations
    if len(slices) == 1:
        for sequence in _cut_sequences(records, slices[0]):
            yield sequence
    else:
        # For multiple slices, concatenate the slice results
        for record in records:
            pieces = (record[s] for s in slices)
            # SeqRecords support addition as concatenation
            yield reduce(lambda x, y: x + y, pieces)

def _update_slices(record, slices):
    n = itertools.count().next
    # Generate a map from indexes in the specified sequence to those in the
    # alignment
    ungap_map = dict((n(), i) for i, base in enumerate(str(record.seq))
                     if base not in GAP_CHARS)
    def update_slice(s):
        """
        Maps a slice relative to ungapped record_id to a slice valid for the
        whole alignment.
        """
        start, end = s.start, s.stop
        if start is not None:
            try:
                start = ungap_map[start]
            except KeyError:
                raise KeyError("""No index {0} in {1}.""".format(
                    start, record.id))
        if end is not None:
            # We need the base in the slice identified by end, not the base
            # at end, otherwise insertions between end-1 and end will be
            # included.
            try:
                end = ungap_map[end - 1] + 1
            except KeyError:
                logging.warn("""No index %d in %s. Keeping columns to end
                    of alignment.""", end, record.id)
                end = None

        return slice(start, end)

    return [update_slice(s) for s in slices]

def cut_sequences_relative(records, slices, record_id):
    """
    Cuts records to slices, indexed by non-gap positions in record_id
    """
    with _record_buffer(records) as r:
        try:
            record = next(i for i in r() if i.id == record_id)
        except StopIteration:
            raise ValueError("Record with id {0} not found.".format(record_id))

        new_slices = _update_slices(record, slices)
        for record in multi_cut_sequences(r(), new_slices):
            yield record

def multi_mask_sequences(records, slices):
    """
    Replace characters sliced by slices with gap characters.
    """
    for record in records:
        record_indices = range(len(record))
        keep_indices = reduce(lambda i, s: i - frozenset(record_indices[s]),
                              slices, frozenset(record_indices))
        seq = ''.join(b if i in keep_indices else '-'
                      for i, b in enumerate(str(record.seq)))
        record.seq = Seq(seq)
        yield record

def mask_sequences_relative(records, slices, record_id):
    with _record_buffer(records) as r:
        try:
            record = next(i for i in r() if i.id == record_id)
        except StopIteration:
            raise ValueError("Record with id {0} not found.".format(record_id))

        new_slices = _update_slices(record, slices)
        for record in multi_mask_sequences(r(), new_slices):
            yield record


def lower_sequences(records):
    """
    Convert sequences to all lowercase.
    """
    logging.info('Applying _lower_sequences generator: '
                 'converting sequences to all lowercase.')
    for record in records:
        yield record.lower()


def upper_sequences(records):
    """
    Convert sequences to all uppercase.
    """
    logging.info('Applying _upper_sequences generator: '
                 'converting sequences to all uppercase.')
    for record in records:
        yield record.upper()


def prune_empty(records):
    """
    Remove any sequences which are entirely gaps ('-')
    """
    for record in records:
        if not all(c == '-' for c in str(record.seq)):
            yield record


def _reverse_annotations(old_record, new_record):
    """
    Copy annotations form old_record to new_record, reversing any
    lists / tuples / strings.
    """
    # Copy the annotations over
    for k, v in old_record.annotations.items():
        # Trim if appropriate
        if isinstance(v, (tuple, list)) and len(v) == len(old_record):
            assert len(v) == len(old_record)
            v = v[::-1]
        new_record.annotations[k] = v

    # Letter annotations must be lists / tuples / strings of the same
    # length as the sequence
    for k, v in old_record.letter_annotations.items():
        assert len(v) == len(old_record)
        new_record.letter_annotations[k] = v[::-1]


def reverse_sequences(records):
    """
    Reverse the order of sites in sequences.
    """
    logging.info('Applying _reverse_sequences generator: '
                 'reversing the order of sites in sequences.')
    for record in records:
        rev_record = SeqRecord(record.seq[::-1], id=record.id,
                               name=record.name,
                               description=record.description)
        # Copy the annotations over
        _reverse_annotations(record, rev_record)

        yield rev_record


def reverse_complement_sequences(records):
    """
    Transform sequences into reverse complements.
    """
    logging.info('Applying _reverse_complement_sequences generator: '
                 'transforming sequences into reverse complements.')
    for record in records:
        rev_record = SeqRecord(record.seq.reverse_complement(),
                               id=record.id, name=record.name,
                               description=record.description)
        # Copy the annotations over
        _reverse_annotations(record, rev_record)

        yield rev_record


def ungap_sequences(records, gap_chars=GAP_CHARS):
    """
    Remove gaps from sequences, given an alignment.
    """
    logging.info('Applying _ungap_sequences generator: '
                 'removing gaps from the alignment.')
    for record in records:
        yield ungap_all(record, gap_chars)

def ungap_all(record, gap_chars=GAP_CHARS):
    record = SeqRecord(Seq(str(record.seq).translate(None, gap_chars)),
            id=record.id, description=record.description)
    return record

def _update_id(record, new_id):
    """
    Update a record id to new_id, also modifying the ID in record.description
    """
    old_id = record.id
    record.id = new_id

    # At least for FASTA, record ID starts the description
    record.description = re.sub('^' + re.escape(old_id), new_id,
            record.description)
    return record


def name_append_suffix(records, suffix):
    """
    Given a set of sequences, append a suffix for each sequence's name.
    """
    logging.info('Applying _name_append_suffix generator: '
                 'Appending suffix ' + suffix + ' to all '
                 'sequence IDs.')
    for record in records:
        new_id = record.id + suffix
        _update_id(record, new_id)
        yield record


def name_insert_prefix(records, prefix):
    """
    Given a set of sequences, insert a prefix for each sequence's name.
    """
    logging.info('Applying _name_insert_prefix generator: '
                 'Inserting prefix ' + prefix + ' for all '
                 'sequence IDs.')
    for record in records:
        new_id = prefix + record.id
        _update_id(record, new_id)
        yield record



def name_include(records, filter_regex):
    """
    Given a set of sequences, filter out any sequences with names
    that do not match the specified regular expression.  Ignore case.
    """
    logging.info('Applying _name_include generator: '
                 'including only IDs matching ' + filter_regex +
                 ' in results.')
    regex = re.compile(filter_regex, re.I)
    for record in records:
        if regex.search(record.id) or regex.search(record.description):
            yield record


def name_exclude(records, filter_regex):
    """
    Given a set of sequences, filter out any sequences with names
    that match the specified regular expression.  Ignore case.
    """
    logging.info('Applying _name_exclude generator: '
                 'excluding IDs matching ' + filter_regex + ' in results.')
    regex = re.compile(filter_regex, re.I)
    for record in records:
        if not regex.search(record.id) and not regex.search(record.description):
            yield record


def name_replace(records, search_regex, replace_pattern):
    """
    Given a set of sequences, replace all occurrences of search_regex
    with replace_pattern. Ignore case.
    """
    regex = re.compile(search_regex, re.I)
    for record in records:
        record.id = regex.sub(replace_pattern, record.id)
        record.description = regex.sub(replace_pattern, record.description)
        yield record


def seq_include(records, filter_regex):
    """
    Filter any sequences who's seq does not match the filter. Ignore case.
    """
    regex = re.compile(filter_regex, re.I)
    for record in records:
        if regex.search(str(record.seq)):
            yield record


def seq_exclude(records, filter_regex):
    """
    Filter any sequences who's seq matches the filter. Ignore case.
    """
    regex = re.compile(filter_regex, re.I)
    for record in records:
        if not regex.search(str(record.seq)):
            yield record


def head(records, head):
    """
    Limit results to the top N records.
    """
    logging.info('Applying _head generator: '
                 'limiting results to top ' + str(head) + ' records.')
    return itertools.islice(records, head)


def tail(records, tail):
    """
    Limit results to the bottom N records.
    """
    with _record_buffer(records) as r:
        record_count = sum(1 for record in r())
        start_index = record_count - tail
        rec_iter = r()
        for record in itertools.islice(rec_iter, start_index, None):
            yield record

# Squeeze-related
def gap_proportion(sequences, gap_chars='-'):
    """
    Generates a list with the proportion of gaps by index in a set of
    sequences.
    """
    aln_len = None
    gaps = []
    for i, sequence in enumerate(sequences):
        if aln_len is None:
            aln_len = len(sequence)
            gaps = [0] * aln_len
        else:
            if not len(sequence) == aln_len:
                raise ValueError(("Unexpected sequence length {0}. Is this "
                                  "an alignment?").format(len(sequence)))

        # Update any gap positions in gap list
        for j, char in enumerate(sequence.seq):
            if char in gap_chars:
                gaps[j] += 1

    sequence_count = float(i + 1)
    gap_props = [i / sequence_count for i in gaps]
    return gap_props


def squeeze(records, gap_threshold=1.0):
    """
    Remove any gaps that are present in the same position across all sequences
    in an alignment.  Takes a second sequence iterator for determining gap
    positions.
    """
    with _record_buffer(records) as r:
        gap_proportions = gap_proportion(r())

        keep_columns = [g < gap_threshold for g in gap_proportions]

        for record in r():
            sequence = str(record.seq)
            # Trim
            squeezed = itertools.compress(sequence, keep_columns)
            yield SeqRecord(Seq(''.join(squeezed)), id=record.id,
                            description=record.description)

def strip_range(records):
    """
    Cut off trailing /<start>-<stop> ranges from IDs.  Ranges must be 1-indexed and
    the stop integer must not be less than the start integer.
    """
    logging.info('Applying _strip_range generator: '
                 'removing /<start>-<stop> ranges from IDs')
    # Split up and be greedy.
    cut_regex = re.compile(r"(?P<id>.*)\/(?P<start>\d+)\-(?P<stop>\d+)")
    for record in records:
        name = record.id
        match = cut_regex.match(str(record.id))
        if match:
            sequence_id = match.group('id')
            start = int(match.group('start'))
            stop = int(match.group('stop'))
            if start > 0 and start <= stop:
                name = sequence_id
        yield SeqRecord(record.seq, id=name,
                        description='')


def transcribe(records, transcribe):
    """
    Perform transcription or back-transcription.
    transcribe must be one of the following:
        dna2rna
        rna2dna
    """
    logging.info('Applying _transcribe generator: '
                 'operation to perform is ' + transcribe + '.')
    for record in records:
        sequence = str(record.seq)
        description = record.description
        name = record.id
        if transcribe == 'dna2rna':
            dna = Seq(sequence, IUPAC.ambiguous_dna)
            rna = dna.transcribe()
            yield SeqRecord(rna, id=name, description=description)
        elif transcribe == 'rna2dna':
            rna = Seq(sequence, IUPAC.ambiguous_rna)
            dna = rna.back_transcribe()
            yield SeqRecord(dna, id=name, description=description)

# Translate-related functions
class CodonWarningTable(object):
    """
    Translation table for codons tht prints a warning when an unknown
    codon is requested, then returns the value passed as missing_char
    """

    def __init__(self, wrapped, missing_char='X'):
        self.wrapped = wrapped
        self.missing_char = missing_char
        self.seen = set()

    def get(self, codon, missing=None):
        try:
            return self.__getitem__(codon)
        except KeyError:
            return missing

    def __getitem__(self, codon):
        if codon == '---':
            return '-'
        elif '-' in codon:
            if codon not in self.seen:
                logging.warn("Unknown Codon: %s", codon)
                self.seen.add(codon)
            return self.missing_char
        else:
            return self.wrapped.__getitem__(codon)

def translate(records, translate):
    """
    Perform translation from generic DNA/RNA to proteins.  Bio.Seq
    does not perform back-translation because the codons would
    more-or-less be arbitrary.  Option to translate only up until
    reaching a stop codon.  translate must be one of the following:
        dna2protein
        dna2proteinstop
        rna2protein
        rna2proteinstop
    """
    logging.info('Applying translation generator: '
                 'operation to perform is ' + translate + '.')

    to_stop = translate.endswith('stop')

    source_type = translate[:3]
    alphabet = {'dna': IUPAC.ambiguous_dna, 'rna': IUPAC.ambiguous_rna}[source_type]

    # Get a translation table
    table = {'dna': CodonTable.ambiguous_dna_by_name["Standard"],
             'rna': CodonTable.ambiguous_rna_by_name["Standard"]}[source_type]

    # Handle ambiguities by replacing ambiguous codons with 'X'
    table.forward_table = CodonWarningTable(table.forward_table)

    for record in records:
        sequence = str(record.seq)
        seq = Seq(sequence, alphabet)
        protein = seq.translate(table, to_stop=to_stop)
        yield SeqRecord(protein, id=record.id, description=record.description)


def max_length_discard(records, max_length):
    """
    Discard any records that are longer than max_length.
    """
    logging.info('Applying _max_length_discard generator: '
                 'discarding records longer than '
                 '.')
    for record in records:
        if len(record) > max_length:
            # Discard
            logging.debug('Discarding long sequence: %s, length=%d',
                record.id, len(record))
        else:
            yield record


def min_length_discard(records, min_length):
    """
    Discard any records that are shorter than min_length.
    """
    logging.info('Applying _min_length_discard generator: '
                 'discarding records shorter than %d.', min_length)
    for record in records:
        if len(record) < min_length:
            logging.debug('Discarding short sequence: %s, length=%d',
                record.id, len(record))
        else:
            yield record


def min_ungap_length_discard(records, min_length):
    """
    Discard any records that are shorter than min_length after removing gaps.
    """
    for record in records:
        if len(ungap_all(record)) >= min_length:
            yield record


def sort_length(source_file, source_file_type, direction=1):
    """
    Sort sequences by length. 1 is ascending (default) and 0 is descending.
    """
    direction_text = 'ascending' if direction == 1 else 'descending'

    logging.info('Indexing sequences by length: %s', direction_text)

    # Adapted from the Biopython tutorial example.

    #Get the lengths and ids, and sort on length
    len_and_ids = sorted((len(rec), rec.id)
                         for rec in SeqIO.parse(source_file, source_file_type))

    if direction == 0:
        ids = reversed([seq_id for (length, seq_id) in len_and_ids])
    else:
        ids = [seq_id for (length, seq_id) in len_and_ids]
    del len_and_ids #free this memory
    record_index = SeqIO.index(source_file.name, source_file_type)
    records = (record_index[seq_id] for seq_id in ids)

    return records


def sort_name(source_file, source_file_type, direction=1):
    """
    Sort sequences by name. 1 is ascending (default) and 0 is descending.
    """
    direction_text = 'ascending' if direction == 1 else 'descending'

    logging.info("Indexing sequences by name: %s", direction_text)

    # Adapted from the Biopython tutorial example.

    #Sort on id
    ids = sorted((rec.id) for rec in SeqIO.parse(source_file,
                                                 source_file_type))

    if direction == 0:
        ids = reversed(ids)
    record_index = SeqIO.index(source_file.name, source_file_type)
    records = (record_index[id] for id in ids)

    return records
