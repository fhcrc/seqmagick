"""
"""
import csv
import os
import os.path
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqIO import FastaIO

import transform
from fileformat import lookup_file_type


class MagickWrap(object):
    """
    A class that wraps functionality present in BioPython.
    """

    def __init__(self, tmp_dir, in_files, out_file=None, alphabet=None,
                 debug=False, verbose=False):
        """
        Constructor
        """
        self.source_files = in_files
        self.tmp_dir = tmp_dir
        self.destination_file = out_file
        self.debug = debug
        self.verbose = verbose

    # Public Methods
    def describe_sequence_files(self, output_format, width):
        """
        Given one more more sequence files, determine if the file is an
        alignment, the maximum sequence length and the total number of
        sequences.  Provides different output formats including tab
        (tab-delimited), csv and align (aligned as if part of a borderless
        table).
        """
        handle = sys.stdout
        if self.destination_file:
            handle = open(self.destination_file, 'w')

        # Create and write out the header row.
        header = ['name', 'alignment', 'min_len', 'max_len', 'avg_len',
                  'num_seqs']
        self._print_file_info(header, output_format=output_format,
                              handle=handle, width=width)
         # Go through all source files passed in, one by one.
        for source_file in self.source_files:
            is_alignment = True
            avg_length = None
            min_length = sys.maxint
            max_length = 0
            sequence_count = 0
            source_file_type = lookup_file_type(os.path.splitext(source_file)[1])

            # Get an iterator and analyze the data.
            for record in SeqIO.parse(source_file, source_file_type):
                # We've found another sequence...
                sequence_count += 1
                sequence_length = len(record)
                if max_length != 0:
                    # If even one sequence is not the same length as the others,
                    # we don't consider this an alignment.
                    if sequence_length != max_length:
                        is_alignment = False

                # Work on determining the length of the longest sequence.
                if sequence_length > max_length:
                    max_length = sequence_length
                if sequence_length < min_length:
                    min_length = sequence_length

                # Average length
                if sequence_count == 1:
                    avg_length = float(sequence_length)
                else:
                    avg_length = avg_length + ((sequence_length - avg_length) /
                                               sequence_count)

            # Handle an empty file:
            if avg_length is None:
                min_length = max_length = avg_length = 0

            self._print_file_info(row=[source_file,
                                       str(is_alignment).upper(),
                                       str(min_length),
                                       str(max_length),
                                       '{0:.2f}'.format(avg_length),
                                       str(sequence_count),
                                       ], output_format=output_format,
                                  handle=handle, width=width)
        if self.destination_file:
            handle.close()

    def create_muscle_alignment(self):
        """
        Use BioPython muscle wrapper to create an alignment.
        """
        muscle_command = MuscleCommandline(input=self.source_files[0], out=self.destination_file)
        if self.debug: print 'DEBUG: muscle command:\n' + str(muscle_command)
        child = subprocess.Popen(str(muscle_command),
                                 stdin=None,
                                 stdout=None,
                                 stderr=None,
                                 shell=(sys.platform!="win32"))
       	return_code = child.wait()
       	return return_code

    def transform(self, cut=False, dash_gap=False, ungap=False, lower=False,
            reverse=False, strict=False, translate=False, upper=False,
            line_wrap=False, first_name_capture=False,
            deduplicate_sequences=False, deduplicate_taxa=False,
            reverse_complement=False, pattern_include=False,
            pattern_exclude=False, squeeze=False, head=False, tail=False,
            sort=False, strip_range=False, transcribe=False, max_length=False,
            min_length=False, name_prefix=False, name_suffix=False,
            input_format=None, output_format=None, prune_empty=False,
            pattern_replace=None, seq_pattern_include=None,
            seq_pattern_exclude=None):
        """
        This method wraps many of the transformation generator functions found
        in this class.
        """

        for source_file in self.source_files:
            # Get just the file name, useful for naming the temporary file.
            file_name = os.path.split(source_file)[1]
            file_ext = os.path.splitext(source_file)[1]
            source_file_type = (input_format or lookup_file_type(file_ext))

            # Specify full path to temporary file for operations that require this.
            # tmp_file will have a seqmagick prefix, i.e. /tmp/seqmagick.a.fasta.
            # If destination_file is part of the magickwrap instance, use that instead
            if self.destination_file is not None:
                destination_file = self.destination_file
            else:
                # Generate a named temporary file
                with tempfile.NamedTemporaryFile(prefix='seqmagick.',
                                                 suffix=file_name,
                                                 delete=False,
                                                 dir=self.tmp_dir) as t:
                    destination_file = t.name

            output_ext = os.path.splitext(destination_file)[1]
            destination_file_type = (output_format or lookup_file_type(output_ext))

            # Get an iterator.
            if sort:
                # Sorted iterator.
                if sort == 'length-asc':
                    records = transform.sort_length(source_file=source_file,
                                                source_file_type=source_file_type,
                                                direction=1)
                elif sort == 'length-desc':
                    records = transform.sort_length(source_file=source_file,
                                                source_file_type=source_file_type,
                                                direction=0)
                elif sort == 'name-asc':
                    records = transform.sort_name(source_file=source_file,
                                              source_file_type=source_file_type,
                                              direction=1)
                elif sort == 'name-desc':
                    records = transform.sort_name(source_file=source_file,
                                              source_file_type=source_file_type,
                                              direction=0)

            else:
                # Unsorted iterator.
                records = SeqIO.parse(source_file, source_file_type)


            #########################################
            # Apply generator functions to iterator.#
            #########################################

            if self.verbose:
                print 'Setting up generator functions for file: ' + source_file

            # Deduplication occurs first, to get a checksum of the
            # original sequence and to store the id field before any
            # transformations occur.

            if max_length:
                records = transform.max_length_discard(records, max_length)

            if min_length:
                records = transform.min_length_discard(records, min_length)

            if deduplicate_sequences:
                records = transform.deduplicate_sequences(records)

            if deduplicate_taxa:
                records = transform.deduplicate_taxa(records)

            if dash_gap:
                records = transform.dashes_cleanup(records)

            if first_name_capture:
                records = transform.first_name_capture(records)
            if upper:
                records = transform.upper_sequences(records)

            if lower:
                records = transform.lower_sequences(records)

            if prune_empty:
                records = transform.prune_empty(records)

            if reverse:
                records = transform.reverse_sequences(records)

            if reverse_complement:
                records = transform.reverse_complement_sequences(records)

            if ungap:
                records = transform.ungap_sequences(records)

            if name_prefix:
                records = transform.name_insert_prefix(records, name_prefix)

            if name_suffix:
                records = transform.name_append_suffix(records, name_suffix)

            if pattern_include:
                records = transform.name_include(records, pattern_include)

            if pattern_exclude:
                records = transform.name_exclude(records, pattern_exclude)

            if pattern_replace:
                search_pattern, replace_pattern = pattern_replace
                records = transform.name_replace(records, search_pattern,
                        replace_pattern)

            if head and tail:
                raise ValueError("Error: head and tail are mutually exclusive "
                        "at the moment.")

            if head:
                records = transform.head(records, head)

            if strip_range:
                records = transform.strip_range(records)

            if tail:
                # To know where to begin including records for tail, we need to count
                # the total number of records, which requires going through the entire
                # file and additional time.
                record_count = sum(1 for record in SeqIO.parse(source_file, source_file_type))
                records = transform.tail(records, tail, record_count)

            if transcribe:
                records = transform.transcribe(records, transcribe)

            if translate:
                records = transform.translate(records, translate)

            if squeeze:
                if self.verbose:
                    print 'Performing squeeze, which requires a new iterator for the first pass.'
                gaps = []
                # Need to iterate an additional time to determine which
                # gaps are share between all sequences in an alignment.
                for record in SeqIO.parse(source_file, source_file_type):
                    # Use numpy to prepopulate a gaps list.
                    if len(gaps) == 0:
                        gaps_length = len(record.seq)
                        #gaps = list(ones( (gaps_length), dtype=int16 ))
                        gaps = [1] * gaps_length
                    gaps = map(transform.gap_check, gaps, list(str(record.seq)))
                records = transform.squeeze(records, gaps)
                if self.verbose: print 'List of gaps to remove for alignment created by squeeze.'
                if self.debug: print 'DEBUG: squeeze gaps list:\n' + str(gaps)

            # cut needs to go after squeeze or the gaps list will no longer be relevent.
            # It is probably best not to use squeeze and cut together in most cases.
            if cut:
                records = transform.cut_sequences(records, start=cut[0], end=cut[1])

           # Only the fasta format is supported, as SeqIO.write does not have a 'wrap' parameter.
            if line_wrap is not None and destination_file_type == 'fasta' and source_file_type == 'fasta':
                if self.verbose: print 'Attempting to write out fasta file with linebreaks set to ' + str(line_wrap) + '.'
                with open(destination_file,"w") as handle:
                    writer = FastaIO.FastaWriter(handle, wrap=line_wrap)
                    writer.write_file(records)
            else:
            # Mogrify requires writing all changes to a temporary file by default,
            # but convert uses a destination file instead if one was specified. Get
            # sequences from an iterator that has generator functions wrapping it.
            # After creation, it is then copied back over the original file if all
            # tasks finish up without an exception being thrown.  This avoids
            # loading the entire sequence file up into memory.
                if self.verbose: print 'Read through iterator and write out transformations to file: ' + destination_file
                SeqIO.write(records, destination_file, destination_file_type)

            # Overwrite original file with temporary file, if necessary.
            if self.destination_file is None:
                if self.verbose: print 'Moving temporary file: ' + destination_file + ' back to file: ' + source_file
                shutil.move(destination_file, source_file)

    def _print_file_info(self, row, output_format, handle, width):
        """
        Write out information that describes a sequence file.
        """
        if 'tab' in output_format:
            handle.write("\t".join(row) + "\n")
        elif 'csv' in output_format:
            writer = csv.writer(handle, delimiter=',', quotechar='"',
                                quoting=csv.QUOTE_NONNUMERIC)
            writer.writerow(row)
        elif 'align' in output_format:
            row = map(lambda s: s.ljust(width)[:width], row)
            handle.write("".join(row) + "\n")


