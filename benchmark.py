#!/usr/bin/env python3
import networkit as nk
import math
import itertools
import time
import csv

import platform
hostname = platform.node()

def power10series(begin, end, steps_per_dec):
    i = 1
    while True:
        value = int(begin * 10**(i / steps_per_dec))
        if (value > end): break
        yield value
        i += 1

TimeoutSeconds = 300
Iterations = 5
Ts = [0.0, 0.5, 0.9]
PLEs = [2.2, 2.5, 3.0]
Ns = list(power10series(1e3, 1e7, 3))
AvgDegs = [10, 100, 1000]

def runNkGen(n, deg, alpha, T):
    gen = nk.generators.HyperbolicGenerator(n=n, k=deg, gamma=alpha, T=T)
    G = gen.generate()
    assert(G.numberOfEdges() == 0) # ensure we do not build the graph datastructure

    sampling_time = gen.getSamplingTimeMS()
    preprocess_time = gen.getPreprocessingTimeMS()
    total_time = gen.getTotalTimeMS()
    num_edges = gen.getNumberOfEdges()

    return [num_edges, sampling_time, preprocess_time, total_time]

with open('nkgen_bench.csv', 'w', newline='') as csvfile:
    fields = ["host", "algo", "iter", "T", "PLE", "n", "deg", "time", "edges", "samplingTime", "preprocessTime", "totalTime"]
    writer = csv.DictWriter(csvfile, fieldnames=fields)
    writer.writeheader()

    for iter, T, ple, deg in itertools.product(range(Iterations), Ts, PLEs, AvgDegs):
        time_ms = 0
        skip = False

        for n in Ns:
            if (3*deg > n):
                continue

            print("Iter: % 2d T: %.1f alpha: %.1f n: % 8d deg: % 4d" % (iter, T, ple, n, deg))

            num_edges = -1
            sampling_time = -1
            preprocess_time = -1
            total_time = -1

            if not skip:
                num_edges, sampling_time, preprocess_time, total_time = runNkGen(n, deg, ple, T)
                skip = (total_time > 1e3 * TimeoutSeconds)
                print("  total time: %d ms %s" % (total_time, "<--- Above Time Threshold; skip larger n" if skip else ""))
            else:
                print("  skipped")

            writer.writerow({'host': hostname, 'algo': 'nkgen', 'iter': iter,
                             'T': T, 'PLE':ple, 'n':n, 'deg':deg,
                             "edges": num_edges,
                             'time':time_ms,"samplingTime": sampling_time, "preprocessTime": preprocess_time, "totalTime": total_time})
            csvfile.flush()

