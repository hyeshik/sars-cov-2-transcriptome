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

import pysam, re

cigar_pat = re.compile('\\d+[A-Z]')

def get_switching_sites(cigar, refstart):
    tokens = cigar_pat.findall(cigar)
    refpos = refstart

    for tok in tokens:
        op = tok[-1]
        length = int(tok[:-1])
        if op == 'N':
            yield refpos, refpos + length
        if op in 'MDN':
            refpos += length

def main(infile, outfile):
    for aln in pysam.AlignmentFile(infile):
        if aln.cigarstring is None or 'N' not in aln.cigarstring:
            continue

        strand = '-' if aln.is_reverse else '+'
        chrom = aln.reference_name

        for i, (sstart, send) in enumerate(get_switching_sites(aln.cigarstring, aln.reference_start)):
            print(chrom, sstart, send,
                  '{}#{}'.format(aln.query_name, i),
                  '0', strand, sep='\t', file=outfile)


if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.stdout)
