#!/usr/bin/env python3
#
# Copyright (c) 2018 Hyeshik Chang
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

from poreplex.fast5_file import Fast5Reader
from sklearn.linear_model import TheilSenRegressor
from concurrent import futures
import h5py
import pandas as pd
import numpy as np

SIGNAL_MEAN_STRIDE = 15
SIGNAL_LENGTH = 30000
SIGNAL_WINDOW_COUNT = SIGNAL_LENGTH // SIGNAL_MEAN_STRIDE
MINIMUM_NONJUMP_POSITIONS = 30

def calculate_scaling_params(events, kmer_mean_levels):
    events = pd.DataFrame(events)
    events['pos'] = events['move'].cumsum()
    jump_positions = events[events['move'] > 1]['pos']
    jump_positions = set(jump_positions - 1) | set(jump_positions)
    nonjump_positions = set(events['pos']) - jump_positions
    if len(nonjump_positions) < MINIMUM_NONJUMP_POSITIONS:
        return

    statelevels = []
    statelevels_jump = []
    for pos, posevents in events.groupby('pos'):
        state = posevents['model_state'].iloc[0]
        if '_' in state:
            continue

        medlevel = posevents['mean'].median()
        if pos in nonjump_positions:
            statelevels.append([medlevel, kmer_mean_levels[state]])
        else:
            statelevels_jump.append([medlevel, kmer_mean_levels[state]])

    statelevels_jump = np.array(statelevels_jump)
    statelevels = np.array(statelevels)
    regr = TheilSenRegressor(random_state=922)
    regr.fit(statelevels[:, 0][:, np.newaxis], statelevels[:, 1])

    return regr.coef_[0], regr.intercept_

def read_raw_signal(handle):
    raw_signal = handle.get_raw_data()
    num_windows = min(SIGNAL_LENGTH, len(raw_signal)) // SIGNAL_MEAN_STRIDE

    meanlevels = raw_signal[:num_windows * SIGNAL_MEAN_STRIDE
                            ].reshape((num_windows, SIGNAL_MEAN_STRIDE)).mean(axis=1)
    if len(meanlevels) < SIGNAL_WINDOW_COUNT:
        meanlevels = np.pad(meanlevels, [SIGNAL_WINDOW_COUNT - len(meanlevels), 0],
                            'constant')
    return meanlevels

def process_read(filename, read_id, kmer_mean_levels):
    try:
        f5 = Fast5Reader(filename, read_id)
    except KeyError:
        return None

    try:
        events = f5.get_basecall()['events']
        scaling_params = calculate_scaling_params(events, kmer_mean_levels)
        if scaling_params is None:
            return

        signal_snippet = read_raw_signal(f5)

        return signal_snippet, np.array(scaling_params)
    finally:
        f5.close()

def runonefile(filename, read_ids, kmer_mean_levels):
    ret = []
    for read_id in read_ids:
        res = process_read(filename, read_id, kmer_mean_levels)
        if res is not None:
            ret.append([read_id, res[1][0], res[1][1]])
    return ret

def process_all(seqsummary_file, kmermodel_file):
    kmermodel = pd.read_table(kmermodel_file)
    kmer_mean_levels = kmermodel.set_index('kmer')['level_mean'].to_dict()

    rrtbl = pd.read_table(seqsummary_file, low_memory=False)

    processes = 36

    with futures.ProcessPoolExecutor(processes) as executor:
        jobs = []
        for filename, row in rrtbl.groupby('filename'):
            read_ids = row['read_id'].tolist()
            job = executor.submit(runonefile, filename, read_ids, kmer_mean_levels)
            jobs.append(job)

        for job in jobs:
            for read_id, scale, shift in job.result():
                print(read_id, scale, shift, sep='\t')

if __name__ == '__main__':
    import sys
    process_all(*sys.argv[1:])

