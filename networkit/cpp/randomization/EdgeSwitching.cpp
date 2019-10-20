// networkit-format

#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/randomization/DegreePreservingShuffle.hpp>
#include <networkit/randomization/EdgeSwitching.hpp>

namespace NetworKit {

void EdgeSwitchingInPlace::run(NetworKit::count numberOfSwitches) {
    if (graph->numberOfEdges() < 2 || numberOfSwitches == 0)
        return;

    if (!hasRun) {
        degreeDistribution = std::discrete_distribution<node>(
            graph->numberOfNodes(), //
            0.0, static_cast<double>(graph->numberOfNodes()),
            [&](double x) { return graph->degree(static_cast<node>(x)); });
    }

    auto &urng = Aux::Random::getURNG();
    Aux::SignalHandler handler;
    while (numberOfSwitches--) {
        handler.assureRunning();
        const auto s1 = degreeDistribution(urng);
        const auto s2 = degreeDistribution(urng);

        // we avoid GraphTools::randomNeighbor to avoid the implicit cost of accessing the
        // Aux::Random::getURNG
        const auto t1 = graph->getIthNeighbor(
            s1, std::uniform_int_distribution<index>{0, graph->degree(s1) - 1}(urng));

        if (s2 == t1 || graph->hasEdge(s2, t1))
            continue;

        const auto t2 = graph->getIthNeighbor(
            s2, std::uniform_int_distribution<index>{0, graph->degree(s2) - 1}(urng));

        if (t1 == t2 || s1 == t2 || graph->hasEdge(s1, t2))
            continue;

        graph->swapEdge(s1, t1, s2, t2);

        ++numberOfSwapsPerformed;
    }

    hasRun = true;
}

EdgeSwitching::EdgeSwitching(const NetworKit::Graph &G, bool doPreprocessing)
    : ownedGraph(doPreprocessing ? DegreePreservingShuffle::shuffleGraph(G) : G) {
    this->graph = &ownedGraph;
}

} // namespace NetworKit
