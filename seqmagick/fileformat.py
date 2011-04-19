'''
'''

# Define mappings in a dictionary with extension : BioPython_file_type.  Defined outside
# of the class so static methods can access it.
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

class FileFormat():
    '''
    A class that maps file extensions to BioPython file formats.
    '''



    @staticmethod
    def lookup_file_type(extension):
        '''
        Convert extension to lower case and look up the corresponding BioPython file type.
        '''
        lower = extension.lower()
        if (lower not in EXTENSION_TO_TYPE):
            raise Exception, "SeqMagick does not know how to handle files with extensions like this: " + extension
        else:
            return EXTENSION_TO_TYPE[lower]


    def __init__(self):
        '''
        Constructor - left for now, but all methods may end up being static.
        '''
        pass



