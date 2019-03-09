#!/usr/bin/env python3
import networkit as nk
import math
import itertools
import time
import csv
import threading
import sys
try:
    import thread
except ImportError:
    import _thread as thread
import platform
hostname = platform.node()

def power10series(begin, end, steps_per_dec):
    i = 1
    while True:
        value = int(begin * 10**(i / steps_per_dec))
        if (value > end): break
        yield value
        i += 1

TimeoutSeconds = 100
Iterations = 5
Ts = [0.0, 0.5, 2.0, 5.0, 10.0]
Alphas = [2.2, 2.5, 3.0]
Ns = list(power10series(1e3, 1e7, 3))
AvgDegs = [10, 100]


with open('nkgen_bench.csv', 'w', newline='') as csvfile:
    fields = ["host", "algo", "iter", "T", "alpha", "n", "deg", "time", "edges", "samplingTime", "preprocessTime", "totalTime"]
    writer = csv.DictWriter(csvfile, fieldnames=fields)
    writer.writeheader()

    for iter, T, alpha, deg in itertools.product(range(Iterations), Ts, Alphas, AvgDegs):
        time_ms = 0

        for n in Ns:
            print("Iter: % 2d T: %.1f alpha: %.1f n: % 8d deg: % 4d" % (iter, T, alpha, n, deg))

            num_edges = -1
            sampling_time = -1
            preprocess_time = -1
            total_time = -1

            if (time_ms >= 0):
                try:
                    timer = threading.Timer(TimeoutSeconds, thread.interrupt_main)
                    timer.start()

                    try:
                        start = time.monotonic()
                        gen = nk.generators.HyperbolicGenerator(n=n, k=deg, gamma=alpha, T=T)
                        G = gen.generate()
                        stop = time.monotonic()
                        assert(G.numberOfEdges() == 0) # ensure we do not build the graph datastructure
                        time_ms = (stop - start) * 1e3

                        sampling_time = gen.getSamplingTimeMS()
                        preprocess_time = gen.getPreprocessingTimeMS()
                        total_time = gen.getTotalTimeMS()
                        num_edges = gen.getNumberOfEdges()

                        print(" ... took: %f ms" % time_ms)

                    finally:
                        timer.cancel()
                except:
                    time_ms = -1
                    print(" ... took too long; cancelled")
            else:
                print(" ... smaller problem timed out: skip")

            writer.writerow({'host': hostname, 'algo': 'nkgen', 'iter': iter,
                             'T': T, 'alpha':alpha, 'n':n, 'deg':deg,
                             "edges": num_edges,
                             'time':time_ms,"samplingTime": sampling_time, "preprocessTime": preprocess_time, "totalTime": total_time})
            csvfile.flush()

