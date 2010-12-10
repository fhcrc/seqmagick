"""
@author: bhodges
"""

import os
import re
import subprocess
import sys
import string

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq, SeqRecord
from Bio.SeqUtils.CheckSum import seguid
from fileformat import FileFormat


class MagickWrap(object):
    """
    A class that wraps functionality present in BioPython.    
    """

# Constructor

    def __init__(self, tmp_dir, in_file, out_file=None, alphabet=None):
        """
        Constructor
        """
        file_name = os.path.split(in_file)[1]

        self.source_file = in_file
        self.source_file_type = FileFormat.lookup_file_type(os.path.splitext(in_file)[1])

        # Specify full path to temporary file for operations that require this.
        # tmp_file will have a seqmagick prefix, i.e. /tmp/seqmagick.a.fasta
        self.tmp_file = os.path.join(tmp_dir, 'seqmagick.' + file_name) 
        self.tmp_file_type = self.source_file_type 

        # Not all operations require a destination file.
        self.destination_file = None
        if out_file is not None:
            self.destination_file = out_file
            self.destination_file_type = FileFormat.lookup_file_type(os.path.splitext(out_file)[1])
            
        # Read in the source file to create a list of SeqRecord objects.
        # This should be fine unless the source file is extremely large.
#        with open(in_file, 'r') as handle:
#            self.source_records = list(SeqIO.parse(handle, self.source_file_type))


# Public Methods


    def convert_format(self):
        """
        Convert input file to a different output format.  This will not work for all formats, 
        e.g. going from fastq to fasta or going from a non-alignment fasta file to phylip would not work.
        """
   
        if self.source_file == self.destination_file:
            raise Exception, "source_file and destination_file cannot be the same file."
        
        if self.destination_file is not None:
           SeqIO.convert(self.source_file, self.source_file_type, self.destination_file, self.destination_file_type) 
        else:
            raise Exception, "An output file was not specified.  Required by the convert action."
        pass


    def ungap_alignment(self):
        """

        """
        pass

    def sort_sequences(self):
        """

        """
        pass


    def wrap_file(self):
        """

        """
        pass


    def cut_sequences(self, start, end):
        """

        """
        pass


    def reverse_sequence_sites(self):
        """

        """
        pass


    def translate_sequences(self):
        """

        """
        pass


    def sequences_to_upper(self):
        """

        """
        pass


    def sequences_to_lower(self):
        """

        """
        pass


    def name_filter(self, filter_regex):
        """
        Given a set of sequences, filter out any sequences with names 
        that do not match the specified regular expression.
        """
        pass


    def sequences_to_lower(self):
        """

        """
        pass


    def is_strict_alphabet(self):
        """

        """
        pass


    def create_muscle_alignment(self):
        """
        Use BioPython muscle wrapper to create an alignment.
        """
        muscle_command = MuscleCommandline(input=self.source_file, out=self.destination_file)
        child = subprocess.Popen(str(muscle_command),
                                 stdin=None,
                                 stdout=None,
                                 stderr=None,
                                 shell=(sys.platform!="win32"))
       	return_code = child.wait()
       	return return_code
        
    def mogrify(self, cut=False, dashgap=False, degap=False, lower=False, 
                reverse=False, strict=False, translate=False, upper=False, wrap=False, 
                first_name_capture=False, deduplicate_sequences=False, deduplicate_taxa=False):
        """
        This method wraps many of the transformation generator functions found 
        in this class.
        """

        # Get an iterator.
        records = SeqIO.parse(self.source_file, self.source_file_type)

        #########################################
        # Apply generator functions to iterator.#
        #########################################

        # Deduplication occurs first, to get a checksum of the 
        # original sequence and to store the id field before any 
        # transformations occur.
     
        if deduplicate_sequences:
            records = self._deduplicate_sequences(records)

        if deduplicate_taxa:
            records = self._deduplicate_taxa(records)

        if dashgap:
            records = self._dashes_cleanup(records)        

        if first_name_capture:
            records = self._first_name_capture(records)


        # Mogrify requires writing all changes to a temporary file first by default, 
        # but use the destination file instead if one was specified. Get
        # sequences from an iterator that has generator functions wrapping it. 
        # After creation, it is then copied back over the original file if all 
        # tasks finish up without an exception being thrown.  This avoids 
        # loading the entire sequence file up into memory.

        if self.destination_file is None:
            SeqIO.write(records, self.tmp_file, self.tmp_file_type)
            # Overwrite original file with temporary file.
            os.rename(self.tmp_file, self.source_file)
        else:
            SeqIO.write(records, self.destination_file, self.destination_file_type)

# Private Methods 


    # Generator Functions
 
    def _dashes_cleanup(self, records):
        """
        Take an alignment and convert any undesirable characters such as ? or ~ to -.
        """
        translation_table = string.maketrans("?~", "--")
        for record in records:
            yield SeqRecord(Seq(str(record.seq).translate(translation_table)), 
                            id=record.id, description=record.description)


    def _deduplicate_sequences(self, records):
        """
        Remove any duplicate records with identical sequences, keep the first 
        instance seen and discard additional occurences.
        """
        checksums = set()
        for record in records:
            checksum = seguid(record.seq)
            if checksum in checksums:
                continue
            checksums.add(checksum)
            yield record

             
    def _deduplicate_taxa(self, records):
        """
        Remove any duplicate records with identical IDs, keep the first 
        instance seen and discard additional occurences.
        """
        taxa = set()
        for record in records:
            # Default to full ID, split if | is found.
            taxid = record.id
            if '|' in record.id:
                taxid = int(record.id.split("|")[0])
            if taxid in taxa:
                continue
            taxa.add(taxid)
            yield record


    def _first_name_capture(self, records):
        """
        Take only the first whitespace-delimited word as the name of the sequence.  
        Essentially removes any extra text from the sequence's description.
        """
        whitespace = re.compile(r'\s+')
        for record in records:
            if whitespace.search(record.description):
                yield SeqRecord(record.seq, id=record.id, 
                                description="")
                                #description=whitespace.split(record.description, 1)[0])
            else: 
                yield record
 
