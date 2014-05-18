"""
Calculate molecular weight, theoretical isoelectric point 
and other physiochemical properties of protein sequences. 
"""
import argparse
import sys

from Bio import SeqIO
from Bio.SeqUtils import ProtParam

from seqmagick import fileformat
from seqmagick import transform

from . import common
def molecular_weight (seq):
    """Calculate MW from a protein sequence.
       
       Forked version of the Bio.SeqUtils.ProtParam version
       to handle 'X' amino acids, and not do weird things with
       water masses.
    """
    # make local dictionary for speed
    MwDict = {}
    # remove a molecule of water from the amino acid weight.
    for i in ProtParam.IUPACData.protein_weights:
        MwDict[i] = ProtParam.IUPACData.protein_weights[i]
    # assign unknown amino acids to an average mass per residue value
    MwDict['X'] = 112.5
    MwDict['B'] = 112.5
    MwDict['U'] = 112.5
    MwDict['J'] = 112.5
    MwDict['O'] = 112.5
    MwDict['Z'] = 112.5
    # gaps and stops have no mass
    MwDict['-'] = 0.0
    MwDict['*'] = 0.0
    MW = 0
    for i in seq:
        MW += MwDict[i]
    return MW

def build_parser(parser):
    parser.add_argument('sequence_file', help="Sequence file",
            type=argparse.FileType('r'))
    parser.add_argument('-o', '--output-file', help="Destination trimmed file",
            type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--source-format', default=None)
    parser.add_argument('-d', '--include-description', action='store_true',
            default=False, help="""Include the sequence description in output
            [default: %(default)s]""")
    parser.add_argument('--sort', dest='sort',
        choices=['length-asc', 'length-desc', 'name-asc', 'name-desc',
                 'mass-asc', 'mass-desc', 'pi-asc', 'pi-desc'],
        help='Sort output based on length, name, molecular weight or pI')
  

def action(arguments):
    common.exit_on_sigpipe()

    # Determine file format for input and output
    source_format = (arguments.source_format or
            fileformat.from_filename(arguments.sequence_file.name))

    with arguments.sequence_file:
        sequences = SeqIO.parse(arguments.sequence_file, source_format)
        
        # sort based on name or length
        sorters = {'length': transform.sort_length,
                   'name': transform.sort_name,}
        directions = {'asc': 1, 'desc': 0}
        if arguments.sort:
            sort_on, direction = arguments.sort.split('-')
            reverse = (direction=='desc')
            if (sort_on == 'length') or (sort_on == 'name'):
              # Sorted iterator
              key, direction = arguments.sort.split('-')
              sequences = sorters[key](source_file=arguments.sequence_file,
                  source_file_type=source_format,
                  direction=directions[direction])
        
        stats = []
        for s in sequences:
            params = ProtParam.ProteinAnalysis(str(s.seq))
          
            stats.append((s, molecular_weight(s.seq), 
                             params.isoelectric_point()))
          
            if arguments.sort and sort_on == 'mass':
              stats = sorted(stats, key=lambda stats: stats[1], reverse=reverse)
            elif arguments.sort and sort_on == 'pi':
              stats = sorted(stats, key=lambda stats: stats[2], reverse=reverse)
        
        if arguments.include_description:
            out = ((s[0].description, s[1], s[2]) for s in stats)
        else:
            out = ((s[0].id, s[1], s[2]) for s in stats)
        
        with arguments.output_file:
            for l in out:
                print >> arguments.output_file, "%s\t%.2f\t%.2f" % (l)
