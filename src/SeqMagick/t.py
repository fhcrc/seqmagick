#! /usr/bin/env python

from magickwrap import MagickWrap

magick = MagickWrap(in_file='examples/aligned.fasta', out_file="test.sth")

magick.convert_format()

#print magick.source_file_type
