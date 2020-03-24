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

import pysam
from itertools import groupby
import re
import sys
import pandas as pd
import glob
from Bio import SeqIO
import gzip

LEADER_TRIMMING_UPPER_LIMIT = 100

cigar_pat = re.compile('\\d+[A-Z]')
input_bam = sys.argv[1]
input_fn = sys.argv[2]
output_fn = sys.argv[3]

def get_sam_reads_by_group(filename):
    get_subject_name = lambda aln: aln.reference_name

    with pysam.AlignmentFile(filename) as alnfile:
        for qname, alns in groupby(alnfile, key=get_subject_name):
            yield qname, list(alns)

def get_alignment_lengths(cigar):
    tokens = cigar_pat.findall(cigar)

    left_clip = int(tokens[0][:-1]) if tokens[0][-1] == 'H' else 0
    right_clip = int(tokens[-1][:-1]) if tokens[-1][-1] == 'H' else 0

    refpos = readpos = 0

    for tok in tokens:
        op = tok[-1]
        length = int(tok[:-1])
        if op in 'MI':
            readpos += length
        if op in 'MDN':
            refpos += length
        if op == 'S':
            raise ValueError('S is not handled.')

    return readpos, refpos, left_clip, right_clip

leader_matches = []
for filename in glob.glob('alignments/VeroInf24h.viral-leader.blast/part-*.bam'):
    #print(filename)
    for qname, alns in get_sam_reads_by_group(filename):
        if len(alns) > 1:
            continue

        for aln in alns:
            readstart = aln.reference_start
            salnwidth, qalnwidth, lclip, rclip = get_alignment_lengths(aln.cigarstring)

            vg_range = lclip, lclip + salnwidth
            rd_range = readstart, readstart + qalnwidth

            leader_matches.append([aln.reference_name, rd_range[0], rd_range[1],
                                   vg_range[0], vg_range[1], lclip, rclip])

leader_matches = \
  pd.DataFrame(leader_matches, columns=['read_id', 'read_start', 'read_end',
                                        'query_start', 'query_end',
                                        'left_clip', 'right_clip'])

high_quality_matches = \
    leader_matches[(leader_matches['read_end'] < LEADER_TRIMMING_UPPER_LIMIT) &
                   (leader_matches['right_clip'] < 1)]

trimming_length = high_quality_matches.set_index('read_id')['read_end'].to_dict()

with gzip.open(output_fn, 'wt') as wrtfile:
    for seq in SeqIO.parse(gzip.open(input_fn, 'rt'), 'fastq'):
        trimlen = trimming_length.get(seq.name)
        if trimlen is not None:
            SeqIO.write([seq[trimlen:]], wrtfile, 'fastq')

