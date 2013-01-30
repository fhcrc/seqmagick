#!/bin/bash

set -e
set -u

NEW_VERSION=$1

sed -i -e 's/__version__ =.*/__version__ = '"'$NEW_VERSION'"'/' seqmagick/__init__.py docs/conf.py
sed -i -e "s/release =.*/release = '$NEW_VERSION'/" docs/conf.py
