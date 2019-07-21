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

    return G, nLeft, nRight

def benchmarkCurveball(G, nTrades, nodeLower, nodeUpper):
    """Returns runtime in seconds"""
    alg = nk.randomization.Curveball(G, False, True)

    iters = 1000

    trades = nk.randomization.CurveballUniformTradeGenerator(nTrades // iters, (nodeLower, nodeUpper))

    start = time.time()

    for i in range(iters):
        alg.run(trades)
        alg.getGraph()

    return 1000 * (time.time() - start)

print("Loading graphs")
datasets = [loadGraph(fn) for fn in glob("netflix1000.txt")]

with open("nk-curveball.csv", "w") as out:
    out.write("algo,iter,nNodes,nEdges,nLeft,nRight,nTrades,timeMS\n")
    for iter in range(30):
        for (G, nLeft, nRight) in datasets:
            for directed in [True, False]:
                nTrades = 10 * nRight
                print("Run with nRight = %d ..." % nRight)
                if directed:
                    duration = benchmarkCurveball(G, nTrades, nLeft, nLeft + nRight)
                else:
                    duration = benchmarkCurveball(G.toUndirected(), nTrades, nLeft, nLeft + nRight)

                print(" ... took %f ms" % duration)

                out.write("nkcb-%s,%d,%d,%d,%d,%d,%d,%f\n" % (
                    "direct" if directed else "undirect",
                    iter, G.numberOfNodes(), G.numberOfEdges(),
                    nLeft, nRight, nTrades, duration
                ))
                out.flush()
