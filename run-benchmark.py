#!/usr/bin/env python3
import networkit as nk
import time
from glob import glob


def loadGraph(filename):
    print(" loading ", filename)
    G = nk.graphio.EdgeListReader(" ", 0, directed=True).read(filename)
    nk.overview(G)

    nLeft = 0
    while G.degreeOut(nLeft) > 0:
        nLeft += 1
    nRight = G.numberOfNodes() - nLeft
    nLeft, nRight

    return G.toUndirected(), nLeft, nRight

def benchmarkCurveball(G, nTrades, nodeLower, nodeUpper):
    """Returns runtime in seconds"""
    alg = nk.randomization.Curveball(G)
    trades = nk.randomization.CurveballUniformTradeGenerator(nTrades, (nodeLower, nodeUpper))

    start = time.time()
    alg.run(trades)

    return 1000 * (time.time() - start)

print("Loading graphs")
datasets = [loadGraph(fn) for fn in glob("netflix100*.txt")]

with open("nk-curveball.csv", "w") as out:
    out.write("iter,nNodes,nEdges,nLeft,nRight,nTrades,timeMS\n")
    for iter in range(30):
        for (G, nLeft, nRight) in datasets:
            nTrades = 10 * nRight
            print("Run with nRight = %d ..." % nRight)
            duration = benchmarkCurveball(G, nTrades, nLeft, nLeft + nRight)
            print(" ... took %f ms" % duration)

            out.write("%d,%d,%d,%d,%d,%d,%f" % (
                iter, G.numberOfNodes(), G.numberOfEdges(),
                nLeft, nRight, nTrades, duration
            ))
            out.flush()
