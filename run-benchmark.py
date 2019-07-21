#!/usr/bin/env python3
import networkit as nk
import time
from glob import glob


def loadGraph(filename):
    print(" loading ", filename)
    G = nk.graphio.NetworkitBinaryReader().read(filename)

    nLeft = sum(int(G.degreeOut(i) > 0) for i in range(G.numberOfNodes()))
    assert(0 < nLeft < G.numberOfNodes())

    return G, nLeft

def benchmarkCurveball(G, nTrades, nodeLower, nodeUpper):
    """Returns runtime in seconds"""
    alg = nk.randomization.Curveball(G, False, True)

    trades = nk.randomization.CurveballUniformTradeGenerator(nTrades, (nodeLower, nodeUpper))

    start = time.time()
    alg.run(trades)
    stop = time.time()

    # Make sure some perturbation took place
    Grand = alg.getGraph()
    pert = nk.distance.EdgeSetDistance(G, Grand).run()
    print(" ... perturbation score: ", pert.perturbationScore())
    assert(G.numberOfEdges() != pert.sizeOfIntersection())

    return 1000 * (stop - start)

print("Loading graphs")
datasets = []
for fn in glob("netflix1000*.nkb"):
    start = time.time()
    G, nLeft = loadGraph(fn)
    mid = time.time()
    datasets.append((G, "l", 0, nLeft))
    print(" ... transpose")
    Gt = G.transpose()
    datasets.append((Gt, "r", nLeft, G.numberOfNodes()))
    end = time.time()
    print(" ... took (%.1f ms and %.1f ms)" % (1000 * (mid - start), 1000 * (end - mid)))


with open("nk-curveball.csv", "w") as out:
    out.write("algo,iter,nNodes,nEdges,activeSide,nActive,nTrades,nRight,timeMS\n")

    for iter in range(10):
        for (G, activeSide, lower, upper) in datasets:
            assert(G.isDirected())
            assert(lower < upper)

            # make [upper, lower) contains all nodes with out-edges
            assert(sum((G.degreeOut(i) for i in range(lower, upper))) == G.numberOfEdges())
            assert(sum((G.degreeIn(i) for i in range(lower, upper))) == 0)

            nRight = G.numberOfNodes() - (lower if lower > 0 else upper)

            for directed in [True, False]:
                nActive = upper - lower
                nTrades = 10 * nActive
                print("Run on [%d, %d) with nRight = %d, directed = %d, side = %s ..." % (lower, upper, nRight, directed, activeSide))
                if directed:
                    duration = benchmarkCurveball(G, nTrades, lower, upper)
                else:
                    duration = benchmarkCurveball(G.toUndirected(), nTrades, lower, upper)

                print(" ... took %f ms" % duration)

                out.write("nkcb-%s,%d,%d,%d,%s,%d,%d,%d,%f\n" % (
                    "direct" if directed else "undirect",
                    iter, G.numberOfNodes(), G.numberOfEdges(),
                    activeSide, nActive, nTrades, nRight, duration
                ))
                out.flush()
