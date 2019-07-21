/*
 * CurveballImpl.h
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
#ifndef RANDOMIZATION_CURVEBALL_IMPL_H
#define RANDOMIZATION_CURVEBALL_IMPL_H

#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>
#include <utility>

#include <tlx/algorithm/random_bipartition_shuffle.hpp>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/randomization/CurveballUniformTradeGenerator.hpp>

#include "AdjacencyList.hpp"
#include "TradeList.hpp"
#include "Helper.hpp"

namespace NetworKit {
namespace CurveballDetails {

using nodepair_vector = std::vector<std::pair<node, node>>;

constexpr node INVALID_NODE = std::numeric_limits<node>::max();
constexpr count LISTROW_END = std::numeric_limits<count>::max();

class CurveballImpl {
    using neighbour_vector = std::vector<node>;

public:
	CurveballImpl(const Graph &G, bool allowSelfLoops = false, bool isBipartite = false) :
        G(G),
        allowSelfLoops(allowSelfLoops),
        isBipartite(isBipartite),
        isDirected(G.isDirected()),
        adjList(G),
        tradeList(G.numberOfNodes()),
        numAffectedEdges(0)
    {}

	void run(const trade_vector &trades) {
        tradeList.initialize(trades);
        initialize();

        Aux::SignalHandler handler;
        auto &urng = Aux::Random::getURNG();
        for (const auto trade : trades) {
            handler.assureRunning();

            // Trade partners u and v
            const node u = std::min(trade.first, trade.second);
            const node v = std::max(trade.first, trade.second);

            process_single_trade(u, v, urng);
        }

        hasRun = true;
    }

    void run(CurveballUniformTradeGenerator& generator) {
        if (!isDirected) {
            auto trades = generator.generate();
            return run(trades);
        }

        initialize();

        Aux::SignalHandler handler;
        auto &urng = Aux::Random::getURNG();
        for (auto i = generator.numberOfTrades(); i; --i) {
            handler.assureRunning();

            // Trade partners u and v
            const auto trade = generator.randomTrade(urng);
            process_single_trade(trade.first, trade.second, urng);
        }

        hasRun = true;
	}

	count getNumberOfAffectedEdges() const {
		assert(hasRun);
		return numAffectedEdges;
	}

	Graph getGraph() const {
	    return adjList.getGraph(G);
	}

	nodepair_vector getEdges() const {
	    return adjList.getEdges();
	}

protected:
	const Graph &G;

	bool hasRun {false};
	bool allowSelfLoops {false};
	bool isBipartite {false};
	bool isDirected {false};

	AdjacencyList adjList;
	TradeList tradeList;
	count maxDegree;
	edgeid numAffectedEdges; // affected half-edges

    neighbour_vector common_neighbours;
    neighbour_vector disjoint_neighbours;

    void initialize() {
        if (!hasRun) {
            maxDegree = G.maxDegree();

            // TODO: Do it in parallel
            G.forEdges([&](node u, node v) {
                update(u, v);
            });

            common_neighbours.reserve(maxDegree);
            disjoint_neighbours.reserve(maxDegree);

        } else if (!isDirected) {
            nodepair_vector edges = adjList.getEdges();

            adjList.restructure();

            for (const auto edge : edges) {
                update(edge.first, edge.second);
            }
        }

    }

    void update(const node a, const node b) {
        if (isDirected) {
            adjList.insertNeighbour(a, b);

        } else {
            const tradeid ta = *(tradeList.getTrades(a));
            const tradeid tb = *(tradeList.getTrades(b));

            if (std::tie(ta, a) <= std::tie(tb, b)) {
                adjList.insertNeighbour(a, b);
            } else {
                adjList.insertNeighbour(b, a);
            }
        }
    }

    void process_single_trade(const node u, const node v, std::mt19937_64& urng) {
        numAffectedEdges += adjList.degreeAt(u);
        numAffectedEdges += adjList.degreeAt(v);

        // Shift the tradeList offset for these two, currently was set to trade_count
        if (!isDirected) {
            tradeList.incrementOffset(u);
            tradeList.incrementOffset(v);
        }

        const auto u_begin = adjList.begin(u);
        const auto v_begin = adjList.begin(v);

        auto u_end = adjList.end(u);
        auto v_end = adjList.end(v);

        // Check whether there exist edges of form (u, v) or (v, u).
        // If this is the case, they are removed from the neighbourhood
        bool edge_between_uv = false;
        bool edge_between_vu = false;

        if (!isBipartite && !allowSelfLoops) {
            auto it = std::find(u_begin, u_end, v);
            if (it != u_end) {
                std::swap(*it, *(--u_end));
                edge_between_uv = true;
            }

            if (isDirected) {
                auto it = std::find(v_begin, v_end, u);
                if (it != v_end) {
                    std::swap(*it, *(--v_end));
                    edge_between_vu = true;
                }
            }
        }

        // Compute common / disjoint neighbors
        CurveballDetails::computeCommonDisjointNeighbour(u_begin, u_end, v_begin, v_end,
            common_neighbours, disjoint_neighbours);

        // Shuffle nodes; [begin, splitter) belong to u, [splitter, end) belong to v
        const auto u_setsize = std::distance(u_begin, u_end) - common_neighbours.size();
        const auto splitter = disjoint_neighbours.cbegin() + u_setsize;
        tlx::random_bipartition_shuffle(disjoint_neighbours.begin(), disjoint_neighbours.end(),
            u_setsize, urng);

        if (isDirected) {
            std::copy(common_neighbours.cbegin(), common_neighbours.cend(), u_begin);
            std::copy(common_neighbours.cbegin(), common_neighbours.cend(), v_begin);

            const auto numCommon = common_neighbours.size();
            std::copy(disjoint_neighbours.cbegin(), splitter, u_begin + numCommon);
            std::copy(splitter, disjoint_neighbours.cend(),   v_begin + numCommon);

        } else {
            // Reset fst/snd row
            adjList.resetRow(u);
            adjList.resetRow(v);

            // Send u's neighbours
            if (edge_between_uv) update(u, v);
            for (size_t i = 0; i < u_setsize; ++i)
                update(u, disjoint_neighbours[i]);

            // Send v's neighbours
            for (size_t i = u_setsize; i < disjoint_neighbours.size(); ++i)
                update(v, disjoint_neighbours[i]);

            // Send common neighbours
            for (const auto common : common_neighbours) {
                update(u, common);
                update(v, common);
            }

        }
    }

};


} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_CURVEBALL_IMPL_H
