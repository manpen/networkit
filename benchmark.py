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

def power10series(begin, end, steps_per_dec):
    i = 1
    while True:
        value = int(begin * 10**(i / steps_per_dec))
        if (value > end): break
        yield value
        i += 1

TimeoutSeconds = 100
Ts = [0.0, 0.5, 2.0, 5.0, 10.0]
Alphas = [2.2, 2.5, 3.0]
Ns = list(power10series(1e3, 1e7, 3))
AvgDegs = [10, 100]


with open('nkgen_bench.csv', 'w', newline='') as csvfile:
    fields = ["T", "alpha", "n", "deg", "time"]
    writer = csv.DictWriter(csvfile, fieldnames=fields)
    writer.writeheader()

    for T, alpha, deg in itertools.product(Ts, Alphas, AvgDegs):
        time_ms = 0

        for n in Ns:
            print("T: %.1f alpha: %.1f n: % 8d deg: % 4d" % (T, alpha, n, deg))

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
                        print(" ... took: %f ms" % time_ms)

                    finally:
                        timer.cancel()
                except:
                    time_ms = -1
                    print(" ... took too long; cancelled")
            else:
                print(" ... smaller problem timed out: skip")

            writer.writerow({'T':T, 'alpha':alpha, 'n':n, 'deg':deg, 'time':time_ms})
            csvfile.flush()

