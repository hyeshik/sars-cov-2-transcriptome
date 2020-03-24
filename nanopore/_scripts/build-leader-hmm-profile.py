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

from collections import defaultdict
import pickle
import pysam

LEADER_REFERENCE_LIMIT = 65

class AlignmentProfile:

    def __init__(self, limit):
        self.limit = limit
        self.emissions = defaultdict(lambda: defaultdict(int))
        self.transitions = defaultdict(lambda: defaultdict(int))

    def sample(self, pairs, readseq):
        state = 'S'
        refcursor = 0
        limit = self.limit

        for readpos, refpos, refbase in pairs:
            if refpos is not None and refpos >= limit:
                self.transitions[state]['E'] += 1
                state = 'E'
                break

            if readpos is None:
                if refbase is None: # splicing
                    self.transitions[state]['E'] += 1
                    state = 'E'
                    break
                else: # small deletion
                    newstate = 'D{}'.format(refpos + 1)
                    self.transitions[state][newstate] += 1
                    state = newstate
                    refcursor = refpos # not needed?
            elif refpos is None: # insertion or softclipping
                newstate = 'I{}'.format(refcursor + 1)
                self.transitions[state][newstate] += 1
                state = newstate
                self.emissions[state][readseq[readpos]] += 1
            else: # Match or mismatch
                newstate = 'M{}'.format(refpos + 1)
                self.transitions[state][newstate] += 1
                state = newstate
                self.emissions[state][readseq[readpos]] += 1
                refcursor = refpos

    def get_emissions(self):
        return {k: dict(v) for k, v in self.emissions.items()}

    def get_transitions(self):
        return {k: dict(v) for k, v in self.transitions.items()}



alnprofile = AlignmentProfile(LEADER_REFERENCE_LIMIT)
alninput = 'alignments/VeroInf24h.viral_genome.sorted.bam'
import sys
progout = sys.stdout

lineno = 0
for aln in pysam.AlignmentFile(alninput, reference_filename='refs/SARS-CoV-2.fa'):
    lineno += 1
    if lineno % 1000 == 0:
        if lineno % 10000 == 0:
            progout.write('({})'.format(lineno))
        progout.write('.')
        progout.flush()

    if aln.flag & 20:
        continue

    try:
        alnpairs = aln.get_aligned_pairs(with_seq=True)
    except TypeError:
        progout.write('x')
    else:
        alnprofile.sample(alnpairs, aln.query_sequence)

print(' Done.')

pickle.dump({
    'emissions': alnprofile.get_emissions(),
    'transitions': alnprofile.get_transitions()},
    open('stats/viral_leader_hmm_{}.pickle'.format(LEADER_REFERENCE_LIMIT), 'wb'))

