"""
@author: bhodges
"""

import os

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from fileformat import FileFormat


class MagickWrap(object):
    """
    A class that wraps functionality present in BioPython.    
    """

# Constructor

    def __init__(self, in_file, out_file=None, alphabet=None):
        """
        Constructor
        """

        self.source_file = in_file
        self.source_file_type = FileFormat.lookup_file_type(os.path.splitext(in_file)[1])

        # Not all operations require a destination file.
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


    def alignment_dashes(self):
        """
        Take an alignment and convert any . or : characters to -.
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

    def align_sequences(self):
        """
        Use BioPython muscle wrapper to create an alignment.
        """
        pass


# Private Methods

    def _read_sequences(self):
        """
        
        """
        pass



        
