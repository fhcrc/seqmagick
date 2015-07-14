#!/bin/bash

subcommands="convert info primer-trim quality-filter extract-ids backtrans-align"

for cmd in $subcommands; do
    seqmagick $cmd --help 2>&1 > ${cmd//-/_}.help
done
