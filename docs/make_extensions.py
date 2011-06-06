#!/usr/bin/env python

from seqmagick import fileformat

extension_map = fileformat.EXTENSION_TO_TYPE
items = extension_map.items()
items.sort()

max_key_length = max((len('Extension'), max(len(k) for k in extension_map.keys())))
max_val_length = max((len('Format'), max(len(v) for v in extension_map.values())))
format_string = '{0:' + str(max_key_length) + 's} {1:' + str(max_val_length) + 's}'

with open('extensions.rst', 'w') as fp:
    def print_row(k, v):
        print >> fp, format_string.format(k, v)
    print_row('=' * max_key_length, '=' * max_val_length)
    print_row('Extension', 'Format')
    print_row('=' * max_key_length, '=' * max_val_length)
    for k, v in items:
        print_row(k, v)
    print_row('=' * max_key_length, '=' * max_val_length)
    print >> fp, ''
