#!/bin/bash 
seqmagick convert --include-from-file selection.txt \
  ../aligned.fasta filtered.fasta
