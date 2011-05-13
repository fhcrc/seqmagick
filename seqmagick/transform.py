"""
Functions to transform / filter sequences
"""
import collections
import itertools
import logging
import re
import string

from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_rna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid


def dashes_cleanup(records):
    """
    Take an alignment and convert any undesirable characters such as ? or ~ to
    -.
    """
    logging.info("Applying _dashes_cleanup: converting any ? or ~ to -.")
    translation_table = string.maketrans("?~", "--")
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


def cut_sequences(records, start, end):
    """
    Cut sequences given a one-based range.  Includes last item.
    """
    start = start - 1
    for record in records:
        cut_record = SeqRecord(record.seq[start:end], id=record.id,
                               name=record.name, description=record.description)

        # Copy the annotations over
        for k, v in record.annotations.items():
            # Trim if appropriate
            if isinstance(v, (tuple, list)) and len(v) == len(record):
                v = v[start:end]
            cut_record.annotations[k] = v

        # Letter annotations must be lists / tuples / strings of the same
        # length as the sequence
        for k, v in record.letter_annotations.items():
            cut_record.letter_annotations[k] = v[start:end]

        yield cut_record


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


def ungap_sequences(records):
    """
    Remove gaps from sequences, given an alignment.
    """
    logging.info('Applying _ungap_sequences generator: '
                 'removing gaps from the alignment.')
    for record in records:
        yield SeqRecord(record.seq.ungap("-"), id=record.id,
                        description=record.description)


def name_append_suffix(records, suffix):
    """
    Given a set of sequences, append a suffix for each sequence's name.
    """
    logging.info('Applying _name_append_suffix generator: '
                 'Appending suffix ' + suffix + ' to all '
                 'sequence IDs.')
    for record in records:
        yield SeqRecord(record.seq, id=record.id+suffix,
                        description=record.description)


def name_insert_prefix(records, prefix):
    """
    Given a set of sequences, insert a prefix for each sequence's name.
    """
    logging.info('Applying _name_insert_prefix generator: '
                 'Inserting prefix ' + prefix + ' for all '
                 'sequence IDs.')
    for record in records:
        yield SeqRecord(record.seq, id=prefix+record.id,
                        description=record.description)



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
        if regex.search(record.id):
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
        if not regex.search(record.id):
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


def tail(records, tail, record_count):
    """
    Limit results to the bottom N records.
    """
    logging.info('Applying _tail generator: '
                 'limiting results to bottom ' + str(tail) + ' records.')
    position = 0
    start = record_count - tail
    for record in records:
        if position < start:
            position += 1
            continue
        else:
            yield record


def squeeze(records, gaps):
    """
    Remove any gaps that are present in the same position across all sequences in an alignment.
    """
    logging.info('Applying _squeeze generator: '
                 'removing any gaps that are present '
                 'in the same position across all sequences in an alignment.')
    sequence_length = len(gaps)
    for record in records:
        sequence = list(str(record.seq))
        if len(sequence) != sequence_length:
            raise ValueError("Unexpected sequence length: {0} != {1} "
                    "Is this an alignment?".format(len(sequence),
                                                   sequence_length))
        squeezed = []
        position = 0
        while (position < sequence_length):
            if bool(gaps[position]) is False:
                squeezed.append(sequence[position])
            position += 1
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
            dna = Seq(sequence, generic_dna)
            rna = dna.transcribe()
            yield SeqRecord(rna, id=name, description=description)
        elif transcribe == 'rna2dna':
            rna = Seq(sequence, generic_rna)
            dna = rna.back_transcribe()
            yield SeqRecord(dna, id=name, description=description)


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
    logging.info('Applying _translation generator: '
                 'operation to perform is ' + translate + '.')
    for record in records:
        sequence = str(record.seq)
        description = record.description
        name = record.id
        if translate in ('dna2protein', 'dna2proteinstop'):
            to_stop = translate == 'dna2proteinstop'
            dna = Seq(sequence, generic_dna)
            protein = dna.translate(to_stop=to_stop)
            yield SeqRecord(protein, id=name, description=description)
        elif translate in ('rna2protein', 'rna2proteinstop'):
            to_stop = translate == 'rna2proteinstop'
            rna = Seq(sequence, generic_rna)
            protein = rna.translate(to_stop=to_stop)
            yield SeqRecord(protein, id=name, description=description)


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


# Begin squeeze-related functions

def _is_gap(character):
    """
    """
    if character == '-':
        return 1
    else:
        return 0


def gap_check(gap, character):
    """
    Build up a gaps list that is used on all sequences
    in an alignment.
    """
    # Skip any characters that have already been found
    if gap == 0:
        return gap
    return int(bool(gap) & bool(_is_gap(character)))

# End squeeze-related functions

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
    record_index = SeqIO.index(source_file, source_file_type)
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
    record_index = SeqIO.index(source_file, source_file_type)
    records = (record_index[id] for id in ids)

    return records
