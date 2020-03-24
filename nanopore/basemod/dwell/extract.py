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

from glob import glob
import sys
sys.path.append('/blaze/hyeshik/nest/tink/gh/poreplex/master')
from poreplex import fast5_file
import h5py
import pandas as pd
import numpy as np
import os

TARGET_SITE = 29015
LEFT_END = TARGET_SITE - 30
RIGHT_END = TARGET_SITE + 31
FAST5_PREFIX = '../..'

def get_migration_time(filename, start, end, corrected_group):
    read_id = filename.split('/')[-1].split('.')[0]

    with fast5_file.Fast5Reader(filename, read_id) as f5:
        h5 = f5.handle
        try:
            tnode = h5['Analyses/'+corrected_group+'/BaseCalled_template']
        except KeyError:
            return None
        alnattrs = dict(tnode['Alignment'].attrs)

        if alnattrs['mapped_start'] > start or alnattrs['mapped_end'] < end:
            return None

        events = pd.DataFrame(tnode['Events'][()])
        events['pos'] = events.index.to_series() + alnattrs['mapped_start']
        events = events[events['pos'].between(start, end)]
        return events.set_index('pos')['length']

def main(catalog):
    for _, catent in catalog.iterrows():
        read_id = catent['read_id']
        filename = os.path.join(FAST5_PREFIX, catent['filename'])

        ret = get_migration_time(filename, LEFT_END, RIGHT_END,
                                 catent['corrected_group'])
        if ret is not None:
            ret = ret.to_frame('dwell').reset_index()
            ret['sample'] = catent['sample']
            ret['read_id'] = read_id
            ret.to_csv(sys.stdout, sep='\t', header=False, index=False)

if __name__ == '__main__':
    catalog_file = pd.read_csv('dwelltime-targets.txt', sep='\t')
    catrange = int(sys.argv[1]), int(sys.argv[2])

    main(catalog_file.iloc[catrange[0]:catrange[1]])
