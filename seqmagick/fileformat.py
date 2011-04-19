"""
Mappings from file extensions to biopython types
"""

# Define mappings in a dictionary with extension : BioPython_file_type.
EXTENSION_TO_TYPE = {'.aln': 'clustal',
                     '.fa': 'fasta',
                     '.faa': 'fasta',
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
                     '.sff': 'sff',
                     '.sth': 'stockholm',
                     '.sto': 'stockholm',
                     }


class UnknownExtensionError(ValueError):
    pass


def lookup_file_type(extension):
    """
    Convert extension to lower case and look up the corresponding BioPython file type.
    """
    try:
        return EXTENSION_TO_TYPE[extension.lower()]
    except KeyError:
        raise UnknownExtensionError("SeqMagick does not know how to handle " +
                "files with extensions like this: " + extension)
