"""
Integration tests, mostly to ensure that basic commands continue working after
commits.

Tests invoke seqmagick.scripts.cli.main, and compare the produced output to the
expected.
"""
import os.path

data_dir = os.path.join(os.path.dirname(__file__), "data")

def data_path(*args):
    return os.path.join(data_dir, *args)
