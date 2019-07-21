#!/usr/bin/env python3
import networkit as nk
import time
from glob import glob


def loadGraph(filename):
    print(" loading ", filename)
    G = nk.graphio.NetworkitBinaryReader().read(filename)

    nLeft = 0
    while G.degreeOut(nLeft) > 0:
        nLeft += 1
    nRight = G.numberOfNodes() - nLeft
    nLeft, nRight

    return G, nLeft, nRight

def benchmarkCurveball(G, nTrades, nodeLower, nodeUpper):
    """Returns runtime in seconds"""
    alg = nk.randomization.Curveball(G, False, True)

    trades = nk.randomization.CurveballUniformTradeGenerator(nTrades, (nodeLower, nodeUpper))

    start = time.time()

    alg.run(trades)

    return 1000 * (time.time() - start)


def computePerturbation(G, nTrades, nodeLower, nodeUpper):
    alg = nk.randomization.Curveball(G, False, True)

    iters = 10
    steps = nTrades // iters

    trades = nk.randomization.CurveballUniformTradeGenerator(steps, (nodeLower, nodeUpper))

    start = time.time()

    for i in range(iters):
        alg.run(trades.generate())
        G_rand = alg.getGraph()

        assert(G_rand.numberOfNodes() == G.numberOfNodes())
        assert(G_rand.numberOfEdges() == G.numberOfEdges())

        pert = nk.distance.EdgeSetDistance(G, G_rand).run()
        print(steps * (i+1), alg.getNumberOfAffectedEdges(), G.numberOfEdges(), pert.sizeOfIntersection(), pert.sizeOfUnion())

    assert(G.numberOfEdges() != pert.sizeOfIntersection())

    return 1000 * (time.time() - start)

print("Loading graphs")
datasets = []
for fn in glob("netflix1000*.nkb"):
    start = time.time()
    G, nLeft, nRight = loadGraph(fn)
    mid = time.time()
    datasets.append((G, "l", 0, nLeft))
    print(" ... transpose")
    datasets.append((G.transpose(), "r", nLeft, nLeft + nRight))
    end = time.time()
    print(" ... took (%.1f ms and %.1f ms)" % (1000 * (mid - start), 1000 * (end - mid)))

with open("nk-curveball.csv", "w") as out:
    out.write("algo,iter,nNodes,nEdges,activeSide,nActive,nTrades,timeMS\n")

    for iter in range(10):
        for (G, activeSide, lower, upper) in datasets:
            assert(G.isDirected())

            for i in range(lower, upper):
                assert(G.degreeIn(i) == 0)

            for directed in [True, False]:
                nActive = upper - lower
                nTrades = 10 * nActive
                print("Run with nRight = %d ..." % nRight)
                if directed:
                    duration = benchmarkCurveball(G, nTrades, 0, nLeft)
                else:
                    duration = benchmarkCurveball(G.toUndirected(), nTrades, lower, upper)

                print(" ... took %f ms" % duration)

                out.write("nkcb-%s,%d,%d,%d,%s,%d,%d,%f\n" % (
                    "direct" if directed else "undirect",
                    iter, G.numberOfNodes(), G.numberOfEdges(),
                    activeSide, nActive, nTrades, duration
                ))
                out.flush()
