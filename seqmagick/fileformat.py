"""
Mappings from file extensions to biopython types
"""
import gzip
import bz2
import os.path

# Define mappings in a dictionary with extension : BioPython_file_type.
EXTENSION_TO_TYPE = {'.aln': 'clustal',
                     '.afa': 'fasta',
                     '.fa': 'fasta',
                     '.faa': 'fasta',
                     '.fas': 'fasta',
                     '.fasta': 'fasta',
                     '.fastq': 'fastq',
                     '.ffn': 'fasta',
                     '.fna': 'fasta',
                     '.frn': 'fasta',
                     '.gb': 'genbank',
                     '.gbk': 'genbank',
                     '.needle': 'emboss',
                     '.phy': 'phylip',
                     '.phylip': 'phylip',
                     '.phyx': 'phylip-relaxed',
                     '.qual': 'qual',
                     '.sff': 'sff',
                     '.sth': 'stockholm',
                     '.sto': 'stockholm',
                     }

COMPRESS_EXT = {'.bz2': bz2.BZ2File, '.gz': gzip.open}


class UnknownExtensionError(ValueError):
    pass


def from_extension(extension):
    """
    Look up the BioPython file type corresponding with input extension.

    Lookup is case insensitive.
    """
    if not extension.startswith('.'):
        raise ValueError("Extensions must begin with a period.")
    try:
        return EXTENSION_TO_TYPE[extension.lower()]
    except KeyError:
        raise UnknownExtensionError("seqmagick does not know how to handle " +
                "files with extensions like this: " + extension)


def from_filename(file_name):
    """
    Lookup the BioPython file type corresponding to an input file name.
    """
    base, extension = os.path.splitext(file_name)
    if extension in COMPRESS_EXT:
        # Compressed file
        extension = os.path.splitext(base)[1]
    return from_extension(extension)
