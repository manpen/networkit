#!/usr/bin/env python3

import networkit as nk
from time import time as timestamp
from glob import glob

for fn_txt in glob("netflix100*.txt"):
    fn_bin = fn_txt.replace(".txt", ".nkb")
    assert(fn_txt != fn_bin)
    print(fn_txt, " -> ", fn_bin)

    start = timestamp()
    G = nk.graphio.EdgeListReader(" ", 0, directed=True).read(fn_txt)
    read_done = timestamp()
    nk.graphio.NetworkitBinaryWriter().write(G, fn_bin)
    write_done = timestamp()
    G1 = nk.graphio.NetworkitBinaryReader().read(fn_bin)
    readbin_done = timestamp()

    assert(G1.numberOfEdges() == G.numberOfEdges())

    print(read_done - start, write_done - read_done, readbin_done - write_done)
