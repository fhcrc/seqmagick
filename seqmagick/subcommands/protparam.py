"""
Extract the sequence IDs from a file
"""
import argparse
import sys

from Bio import SeqIO
from Bio.SeqUtils import ProtParam

from seqmagick import fileformat

from . import common

def build_parser(parser):
    parser.add_argument('sequence_file', help="Sequence file",
            type=argparse.FileType('r'))
    parser.add_argument('-o', '--output-file', help="Destination trimmed file",
            type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--source-format', default=None)
    parser.add_argument('-d', '--include-description', action='store_true',
            default=False, help="""Include the sequence description in output
            [default: %(default)s]""")
    parser.add_argument('-s','--sort', default=False, 
                        help="Sort output based on molecuar weight",
                        action='store_true')
  

def action(arguments):
    common.exit_on_sigpipe()

    # Determine file format for input and output
    source_format = (arguments.source_format or
            fileformat.from_filename(arguments.sequence_file.name))

    with arguments.sequence_file:
        sequences = SeqIO.parse(arguments.sequence_file, source_format)
        stats = []
        for s in sequences:
          params = ProtParam.ProteinAnalysis(str(s.seq))
          
          stats.append((s, params.molecular_weight()))
        if arguments.sort:
          stats = sorted(stats, key=lambda stats: stats[1])
        
        if arguments.include_description:
            out = ((s[0].description, s[1]) for s in stats)
        else:
            out = ((s[0].id, s[1]) for s in stats)
        
        with arguments.output_file:
            for l in out:
                print >> arguments.output_file, "%s\t%.2f" % (l)
