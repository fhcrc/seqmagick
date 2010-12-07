#! /usr/bin/env python

import sys
from Bio import AlignIO
from Bio import SeqIO

ALIGNMENT_FILE = 'alignment.fasta'
UNGAPPED_FILE = 'ungapped.fasta'

# Ungap all sequences present in an alignment.
def ungap_sequences():
    alignment_file = open(ALIGNMENT_FILE, 'r')
    alignment = AlignIO.read(alignment_file, "fasta")
    alignment_file.close()
    ungapped = []
    for record in alignment:
        record.seq = record.seq.ungap("-")
        ungapped.append(record)
    ungapped_file = open(UNGAPPED_FILE, 'w') 
    SeqIO.write(ungapped, ungapped_file, "fasta")
    ungapped_file.close()

def main():
    ungap_sequences()


if __name__ == '__main__':
    sys.exit(main())

