#! /usr/bin/env python
# --squeeze prototype


def is_gap(c):
    if c == '-':
        return 1
    else:
        return 0


def gap_check(i, c):
    if i == 0:
        return i
    i = int(bool(i) & bool(is_gap(c)))
    return i
    

def squeeze(gaps, sequence):
    squeezed = []
    sequence_length = len(gaps) 
    position = 0
    while (position < sequence_length):
        if bool(gaps[position]) is False:
            squeezed.append(sequence[position])
        position += 1
    return ''.join(squeezed)

gaps = [1, 1, 1, 1, 1]
s1 = ['A', '-', 'A', '-', 'A']
s2 = ['A', 'A', 'A', '-', 'A']
s3 = ['A', 'A', 'A', '-', '-']

for s in [s1, s2, s3]:
    gaps = map(gap_check, gaps, s)

# would be nice to check to make sure all sequences are the same length.

#gaps = map(gap_check, gaps, s1)

#print gaps

for s in [s1, s2, s3]:
    print squeeze(gaps, s)

