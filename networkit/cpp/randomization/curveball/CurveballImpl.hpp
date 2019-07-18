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
	CurveballImpl(const Graph &G) :
        G(G),
        numNodes(G.numberOfNodes()),
        adjList(G),
        tradeList(G.numberOfNodes()),
        numAffectedEdges(0)
    {
        hasRun = false;
        assert(G.checkConsistency());
        assert(G.numberOfSelfLoops() == 0);
        assert(numNodes > 0);
    }

    template <bool Directed = false>
	void run(const trade_vector &trades) {
        if (!hasRun)
            loadFromGraph(trades);
        else
            restructureGraph(trades);

        neighbour_vector common_neighbours;
        neighbour_vector disjoint_neighbours;

        common_neighbours.reserve(maxDegree);
        disjoint_neighbours.reserve(maxDegree);

        Aux::SignalHandler handler;

        auto &urng = Aux::Random::getURNG();

        constexpr bool allowSelfLoops = false;

        for (const auto trade : trades) {
            handler.assureRunning();

            // Trade partners u and v
            const node u = std::min(trade.first, trade.second);
            const node v = std::max(trade.first, trade.second);

            numAffectedEdges += adjList.degreeAt(u);
            numAffectedEdges += adjList.degreeAt(v);

            // Shift the tradeList offset for these two, currently was set to trade_count
            tradeList.incrementOffset(u);
            tradeList.incrementOffset(v);

            const auto u_begin = adjList.begin(u);
            const auto v_begin = adjList.begin(v);

            auto u_end = adjList.end(u);
            auto v_end = adjList.end(v);

            // Check whether there exist edges of form (u, v) or (v, u).
            // If this is the case, they are removed from the neighbourhood
            bool edge_between_uv = false;
            bool edge_between_vu = false;
            if (!Directed || !allowSelfLoops) {
                auto it = std::find(u_begin, u_end, v);
                if (it != u_end) {
                    *it = *(--u_end);
                    edge_between_uv = true;
                }

                if (Directed) {
                    auto it = std::find(v_begin, v_end, u);
                    if (it != v_end) {
                        *it = *(--v_end);
                        edge_between_vu = true;
                    }
                }
            }

            // Compute common / disjoint neighbors
            CurveballDetails::computeCommonDisjointNeighbour(u_begin, u_end, v_begin, v_end, common_neighbours, disjoint_neighbours);

            // Reset fst/snd row
            adjList.resetRow(u);
            adjList.resetRow(v);


            // Shuffle nodes; [begin, splitter) belong to u, [splitter, end) belong to v
            const auto u_setsize = std::distance(u_begin, u_end) - common_neighbours.size();
            const auto splitter = disjoint_neighbours.cbegin() + u_setsize;
            tlx::random_bipartition_shuffle(disjoint_neighbours.begin(), disjoint_neighbours.end(), u_setsize, urng);

            // Send u's neighbours
            if (edge_between_uv) update(u, v);
            for(auto it = disjoint_neighbours.cbegin(); it != splitter; ++it)
                update(u, *it);

            // Send v's neighbours
            if (edge_between_vu) update(v, u);
            for(auto it = splitter; it != disjoint_neighbours.cend(); ++it)
                update(v, *it);

            // Send common neighbours
            for (const auto common : common_neighbours) {
                update(u, common);
                update(v, common);
            }
        }

        hasRun = true;
    }

	count getNumberOfAffectedEdges() const {
		assert(hasRun);
		return numAffectedEdges;
	}

	Graph getGraph() const {
	    return adjList.getGraph();
	}

	nodepair_vector getEdges() const {
	    return adjList.getEdges();
	}

protected:
	const Graph &G;
	const node numNodes;

	bool hasRun;
	AdjacencyList adjList;
	TradeList tradeList;
	count maxDegree;
	edgeid numAffectedEdges; // affected half-edges

	void loadFromGraph(const trade_vector &trades) {
        maxDegree = G.maxDegree();

        tradeList.initialize(trades);

        // Insert to adjacency list, directed according trades
        // TODO: make parallel
        G.forEdges([&](node u, node v) {
            update(u, v);
        });
    }

	void restructureGraph(const trade_vector &trades) {
        // TODO: replace by forEdge
        nodepair_vector edges = adjList.getEdges();

        adjList.restructure();
        tradeList.initialize(trades);

        for (const auto edge : edges) {
            update(edge.first, edge.second);
        }
    }

	void update(const node a, const node b) {
		const tradeid ta = *(tradeList.getTrades(a));
		const tradeid tb = *(tradeList.getTrades(b));

		if (std::tie(ta, a) <= std::tie(tb, b)) {
			adjList.insertNeighbour(a, b);
		} else {
			adjList.insertNeighbour(b, a);
		}
	}
};


} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_CURVEBALL_IMPL_H
