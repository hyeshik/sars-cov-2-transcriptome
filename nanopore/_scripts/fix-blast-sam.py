#!/usr/bin/env python3
#
# Copyright (c) 2020 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

import sys
import gzip
from Bio import SeqIO

SAM_INPUT = sys.argv[1]
QUERY_INPUT = sys.argv[2]
SUBJECT_INPUT = sys.argv[3]

query_ids = []
for seq in SeqIO.parse(QUERY_INPUT, 'fasta'):
    query_ids.append((seq.name, len(seq.seq)))

subject_ids = []
for seq in SeqIO.parse(SUBJECT_INPUT, 'fasta'):
    subject_ids.append((seq.name, len(seq.seq)))

for line in gzip.open(SAM_INPUT, 'rt'):
    fields = line.split('\t')

    if line.startswith('@SQ'):
        qno = int(fields[1][3:].replace('Query_', ''))
        seqid, seqlen = query_ids[qno-1]
        if seqlen != int(fields[2][3:]):
            print('Query_{} != {}'.format(qno, seqid), file=sys.stderr)
        fields[1] = 'SN:' + seqid
    elif line.startswith('@'):
        pass
    else:
        subject_no = int(fields[0].replace('Subject_', '')) - 1
        fields[0] = subject_ids[subject_no][0]
        qno = int(fields[2].replace('Query_', ''))
        fields[2] = query_ids[qno - 1][0]

    print(*fields, sep='\t', end='')


