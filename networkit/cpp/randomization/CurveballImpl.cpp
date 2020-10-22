/*
 * CurveballImpl.cpp
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
// networkit-format
#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

#include <range/v3/all.hpp>

#include "CurveballImpl.hpp"
#include <tlx/algorithm/random_bipartition_shuffle.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>

#include "trades/PermutationTrade.hpp"

namespace NetworKit {
namespace CurveballDetails {

using degree_vector = std::vector<count>;
using trade_vector = std::vector<trade_descriptor>;
using neighbour_vector = std::vector<node>;
using node_vector = std::vector<node>;
using nodepair_vector = std::vector<std::pair<node, node>>;

using neighbour_vector = std::vector<node>;
using degree_vector = std::vector<count>;
using degree_it = std::vector<count>::const_iterator;
using pos_vector = std::vector<edgeid>;
using neighbour_it = neighbour_vector::iterator;
using cneighbour_it = neighbour_vector::const_iterator;
using nodepair_vector = std::vector<std::pair<node, node>>;

/**
 * @brief Initialize method (when constructor can't be used)
 *
 */
void CurveballAdjacencyList::initialize(const degree_vector &degrees, const edgeid degree_count) {
    neighbours.resize(degree_count + degrees.size() + 1);
    offsets.resize(degrees.size());
    begins.resize(degrees.size() + 1);
    degreeCount = degree_count;

    count sum = 0;
    node node_id = 0;
    for (const count node_degree : degrees) {
        begins[node_id] = sum;

        // assert(node_degree > 0);

        sum += node_degree;
        neighbours[sum] = LISTROW_END;

        // shift after Sentinel
        sum += 1;
        node_id++;
    }
    neighbours[sum] = LISTROW_END;
    begins[degrees.size()] = sum;

    assert(sum == degreeCount + degrees.size());
    assert(node_id == degrees.size());

    return;
}

void CurveballAdjacencyList::restructure() {
    std::fill(offsets.begin(), offsets.end(), 0);
    return;
}

/**
 * @brief Constructor
 * @param degree_vector Pointer to a vector with node degrees
 * @param degree_count Sum of all degrees in degree_vector
 *
 * We add to each adjacency list entry a delimiter to mark the end
 */
CurveballAdjacencyList::CurveballAdjacencyList(const degree_vector &degrees,
                                               const edgeid degree_count)
    : neighbours(degree_count + degrees.size() + 1), offsets(degrees.size()),
      begins(degrees.size() + 1), degreeCount(degree_count) {
    count sum = 0;
    node node_id = 0;
    for (const count node_degree : degrees) {
        begins[node_id] = sum;

        // no isolated nodes allowed
        assert(node_degree > 0);

        sum += node_degree;
        neighbours[sum] = LISTROW_END;

        sum += 1;
        node_id++;
    }
    neighbours[sum] = LISTROW_END;
    begins[degrees.size()] = sum;

    assert(sum == degree_count + degrees.size());
    assert(node_id == degrees.size());
}

nodepair_vector CurveballAdjacencyList::getEdges() const {
    nodepair_vector edges;
    edges.reserve(degreeCount);
    for (node nodeid = 0; nodeid < static_cast<node>(offsets.size()); nodeid++) {
        for (auto it = cbegin(nodeid); it != cend(nodeid); it++) {
            edges.push_back(std::make_pair(nodeid, *it));
        }
    }

    return edges;
}

///////////////////////////////////////////////////////////////////////////////

CurveballMaterialization::CurveballMaterialization(const CurveballAdjacencyList &adj_list)
    : adjacencyList(adj_list) {}

Graph CurveballMaterialization::toGraph(bool parallel) {
    Graph G(adjacencyList.numberOfNodes(), false, false);

    if (parallel)
        toGraphParallel(G);
    else
        toGraphSequential(G);

    return G;
}

void CurveballMaterialization::toGraphParallel(Graph &G) {
    // 1) insertion of first half is threadsafe
    // 2) insertion of second half is not threadsafe, for now done sequentially

    const node numNodes = adjacencyList.numberOfNodes();

    std::vector<count> missingEdgesCounts(numNodes, 0);

    std::vector<std::vector<edgeweight>> new_edgeWeights(numNodes);
    std::vector<std::vector<node>> new_outEdges(numNodes);
    std::vector<count> outDeg(numNodes);

    // Add first half of edges and count missing edges for each node
    G.parallelForNodes([&](node nodeid) {
        const count degree =
            static_cast<count>(adjacencyList.cend(nodeid) - adjacencyList.cbegin(nodeid));
        outDeg[nodeid] = degree;
        missingEdgesCounts[nodeid] = adjacencyList.degreeAt(nodeid) - degree;
        new_outEdges[nodeid].reserve(degree);
        new_edgeWeights[nodeid].resize(degree, 1);
        for (auto it = adjacencyList.cbegin(nodeid); it != adjacencyList.cend(nodeid); it++) {
            new_outEdges[nodeid].push_back(*it);
        }
    });

    G.outEdges.swap(new_outEdges);

    // Reserve the space
    G.parallelForNodes([&](node v) { G.outEdges[v].reserve(outDeg[v] + missingEdgesCounts[v]); });

    // Second half of the edges
    G.forNodes([&](node v) {
        for (count neighbor_id = 0; neighbor_id < outDeg[v]; neighbor_id++) {
            const node u = G.outEdges[v][neighbor_id];
            G.outEdges[u].push_back(v);
        }
    });

    // TODO: is the networkit adjacency list even sorted for the neighbours? if
    // not omit this
    // Sort neighbours
    G.parallelForNodes([&](node v) { std::sort(G.outEdges[v].begin(), G.outEdges[v].end()); });

    // Set number of self-loops
    G.storedNumberOfSelfLoops = 0;

    // Set numberOfEdges
    G.m = adjacencyList.numberOfEdges() / 2;

    // Shrink to fit
    G.shrinkToFit();
}

void CurveballMaterialization::toGraphSequential(Graph &G) {
    const node numNodes = adjacencyList.numberOfNodes();

    // Analogue to "toGraphSequential" of GraphBuilder
    std::vector<count> missingEdgesCounts;
    missingEdgesCounts.reserve(numNodes);

    std::vector<std::vector<edgeweight>> new_edgeWeights(numNodes);
    std::vector<std::vector<node>> new_outEdges(numNodes);
    std::vector<count> outDeg(numNodes);

    // Add first half of edges and count missing edges for each node
    G.forNodes([&](node nodeid) {
        const count degree =
            static_cast<count>(adjacencyList.cend(nodeid) - adjacencyList.cbegin(nodeid));
        outDeg[nodeid] = degree;
        missingEdgesCounts.push_back(adjacencyList.degreeAt(nodeid) - degree);
        new_outEdges[nodeid].reserve(degree);
        new_edgeWeights[nodeid].resize(degree, 1);
        for (auto it = adjacencyList.cbegin(nodeid); it != adjacencyList.cend(nodeid); it++) {
            new_outEdges[nodeid].push_back(*it);
        }
    });

    G.outEdges.swap(new_outEdges);

    // Reserve the space
    G.forNodes([&](node v) { G.outEdges[v].reserve(outDeg[v] + missingEdgesCounts[v]); });

    // Second half of the edges
    G.forNodes([&](node v) {
        for (count neighbor_id = 0; neighbor_id < outDeg[v]; neighbor_id++) {
            const node u = G.outEdges[v][neighbor_id];
            G.outEdges[u].push_back(v);
        }
    });

    // TODO: is the networkit adjacency list even sorted for the neighbours? if
    // not omit this
    // Sort neighbours
    G.forNodes([&](node v) { std::sort(G.outEdges[v].begin(), G.outEdges[v].end()); });

    // Set number of self-loops
    G.storedNumberOfSelfLoops = 0;

    // Set numberOfEdges
    G.m = adjacencyList.numberOfEdges() / 2;

    // Shrink to fit
    G.shrinkToFit();
}

///////////////////////////////////////////////////////////////////////////////
TradeList::TradeList(const node num_nodes) : numNodes(num_nodes) {}

void TradeList::initialize(const trade_vector &trades) {
    tradeList.clear();
    tradeList.resize(2 * trades.size() + numNodes);
    offsets.clear();
    offsets.resize(numNodes + 1);

    assert(numNodes > 0);
    assert(!trades.empty());

    std::vector<tradeid> trade_count(numNodes);

    // Push occurrences
    for (const auto &trade : trades) {
        assert(trade.first < numNodes);
        assert(trade.second < numNodes);

        trade_count[trade.first]++;
        trade_count[trade.second]++;
    }

    // add missing +1 for sentinel
    trade_count[0]++;
    std::partial_sum(trade_count.cbegin(), trade_count.cend(), offsets.begin() + 1,
                     [&](const tradeid a, const tradeid b) { return a + b + 1; });
    // add dummy
    offsets[numNodes] = 2 * trades.size() + numNodes - 1;

    // set sentinels
    for (node node = 1; node < numNodes; node++) {
        tradeList[offsets[node] - 1] = TRADELIST_END;
    }
    // set last entry as sentinel
    tradeList.back() = TRADELIST_END;

    std::fill(trade_count.begin(), trade_count.end(), 0);
    {
        tradeid trade_id = 0;
        for (const auto &trade : trades) {
            auto updateNode = [&](const node u) {
                const node pos = offsets[u] + trade_count[u];
                tradeList[pos] = trade_id;
                trade_count[u]++;
            };

            updateNode(trade.first);
            updateNode(trade.second);
            trade_id++;
        }
    }
}

TradeList::TradeList(const trade_vector &trades, const node num_nodes)
    : tradeList(2 * trades.size() + num_nodes), offsets(num_nodes + 1), numNodes(num_nodes) {
    // Manuel: see above

    assert(num_nodes > 0);
    assert(!trades.empty());

    std::vector<tradeid> trade_count(num_nodes);

    // Push occurences
    for (const auto &trade : trades) {
        assert(trade.first < num_nodes);
        assert(trade.second < num_nodes);

        trade_count[trade.first]++;
        trade_count[trade.second]++;
    }

    // add missing +1 for sentinel
    trade_count[0]++;
    std::partial_sum(trade_count.cbegin(), trade_count.cend(), offsets.begin() + 1,
                     [&](const tradeid a, const tradeid b) { return a + b + 1; });
    // add dummy
    offsets[num_nodes] = 2 * trades.size() + num_nodes - 1;

    // set sentinels
    for (node node = 1; node < numNodes; node++) {
        tradeList[offsets[node] - 1] = TRADELIST_END;
    }
    // set last entry as sentinel
    tradeList.back() = TRADELIST_END;

    std::fill(trade_count.begin(), trade_count.end(), 0);
    {
        tradeid trade_id = 0;
        for (const auto &trade : trades) {
            auto updateNode = [&](const node u) {
                const node pos = offsets[u] + trade_count[u];
                tradeList[pos] = trade_id;
                trade_count[u]++;
            };

            updateNode(trade.first);
            updateNode(trade.second);
            trade_id++;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

CurveballIM::CurveballIM(const Graph &G)
    : G(&G), numNodes(G.numberOfNodes()), tradeList(G.numberOfNodes()), numAffectedEdges(0) {
    hasRun = false;
    assert(G.checkConsistency());
    assert(G.numberOfSelfLoops() == 0);
    assert(numNodes > 0);
}

void CurveballIM::loadFromGraph(const std::vector<std::pair<node, node>> &trades) {
    // Compute degree sequence
    degree_vector degrees;
    degrees.reserve(numNodes);
    edgeid degree_sum = 0;
    G->forNodes([&](node v) {
        degrees.push_back(G->degree(v));
        degree_sum += G->degree(v);
    });

    maxDegree = *(std::max_element(degrees.cbegin(), degrees.cend()));

    adjList.initialize(degrees, degree_sum);
    tradeList.initialize(trades);

    // Insert to adjacency list, directed according trades
    G->forEdges([&](node u, node v) { update(u, v); });
    return;
}

void CurveballIM::restructureGraph(const std::vector<std::pair<node, node>> &trades) {
    nodepair_vector edges = adjList.getEdges();

    adjList.restructure();
    tradeList.initialize(trades);

    for (const auto &edge : edges) {
        update(edge.first, edge.second);
    }

    return;
}

void CurveballIM::run(const trade_vector &trades) {
    if (!hasRun)
        loadFromGraph(trades);
    else
        restructureGraph(trades);

    count trade_count = 0;
    neighbour_vector common_neighbours;
    neighbour_vector disjoint_neighbours;

    common_neighbours.reserve(maxDegree);
    disjoint_neighbours.reserve(maxDegree);

    Aux::SignalHandler handler;

    auto &urng = Aux::Random::getURNG();

    Trades::RandomPerumation trader(maxDegree);

    for (const auto &trade : trades) {
        handler.assureRunning();

        // Trade partners u and v
        const node u = trade.first;
        const node v = trade.second;

        numAffectedEdges += adjList.degreeAt(u) + adjList.degreeAt(v);

        // Shift the tradeList offset for these two, currently was set to
        // trade_count
        tradeList.incrementOffset(u);
        tradeList.incrementOffset(v);

        const tradeid tu = *(tradeList.getTrades(u));
        const tradeid tv = *(tradeList.getTrades(v));

        ranges::sort(adjList.range(u));
        ranges::sort(adjList.range(v));

        trader.run(u, adjList.range(u), v, adjList.range(v), urng);

        adjList.resetRow(u);
        adjList.resetRow(v);

        for (auto n : trader.neighborsOfU()) {
            assert(n != u && n != v);
            update(u, tu, n);
        }

        for (auto n : trader.neighborsOfV()) {
            assert(n != u && n != v);
            update(v, tv, n);
        }

        for (auto n : trader.neighborsInCommon()) {
            assert(n != u && n != v);
            update(u, tu, n);
            update(v, tv, n);
        }

        if (trader.hasSharedEdge())
            update(u, tu, v, tv);

        trade_count++;
    }

    hasRun = true;

    return;
}

Graph CurveballIM::getGraph(bool parallel) const {
    CurveballMaterialization gb(adjList);

    return gb.toGraph(parallel);
}

} // namespace CurveballDetails
} // namespace NetworKit
