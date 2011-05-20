"""
A collection of functions to apply
"""
import hashlib

def no_gaps(records):
    for record in records:
        if not '-' in str(record.seq):
            yield record

def hash_starts_numeric(records):
    """
    Very useful function that only accepts records with a numeric start to
    their sha-1 hash.
    """
    for record in records:
        seq_hash = hashlib.sha1(str(record.seq)).hexdigest()
        if seq_hash[0].isdigit():
            yield record
