#!/usr/bin/env python3

import sys

if __name__ == '__main__':
    from seqmagick.scripts import cli
    sys.exit(cli.main(sys.argv[1:]))
