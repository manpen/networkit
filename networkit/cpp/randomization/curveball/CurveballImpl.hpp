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
	CurveballImpl(const Graph &G, bool allowSelfLoops = false, bool isBipartite = false) :
        G(G),
        numNodes(G.numberOfNodes()),
        allowSelfLoops(allowSelfLoops),
        isBipartite(isBipartite),
        adjList(G),
        tradeList(G.numberOfNodes()),
        numAffectedEdges(0)
    {}

	void run(const trade_vector &trades) {
        if (G.isDirected()) {
            process<true>(trades);
        } else {
            process<false>(trades);
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
	const node numNodes;

	bool hasRun {false};
	bool allowSelfLoops {false};
	bool isBipartite {false};

	AdjacencyList adjList;
	TradeList tradeList;
	count maxDegree;
	edgeid numAffectedEdges; // affected half-edges

    template <bool Directed = false>
    void process(const trade_vector& trades) {
        if (!hasRun)
            loadFromGraph<Directed>(trades);
        else
            restructureGraph<Directed>(trades);

        neighbour_vector common_neighbours;
        neighbour_vector disjoint_neighbours;

        common_neighbours.reserve(maxDegree);
        disjoint_neighbours.reserve(maxDegree);

        Aux::SignalHandler handler;

        auto &urng = Aux::Random::getURNG();

        constexpr bool allowSelfLoops = true;

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

            if (!isBipartite && (!Directed || !allowSelfLoops)) {
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

            if (Directed) {
                // Send u's neighbours
                if (edge_between_uv) update<Directed>(u, v);
                for (auto it = disjoint_neighbours.cbegin(); it != splitter; ++it)
                    update<Directed>(u, *it);

                // Send v's neighbours
                if (edge_between_vu) update<Directed>(v, u);
                for (auto it = splitter; it != disjoint_neighbours.cend(); ++it)
                    update<Directed>(v, *it);


                // Send common neighbours
                for (const auto common : common_neighbours) {
                    update<Directed>(u, common);
                    update<Directed>(v, common);
                }

            } else {
                auto prefetched_update = [&] (node u, node* begin, size_t size) {
                    const auto prefetch_len = std::min<size_t>(16, size) / 2;

                    for(size_t i = 0; i < prefetch_len; ++i) {
                        tradeList.prefetchTrades1(begin[i]);
                    }

                    for(size_t i = prefetch_len; i < 2 * prefetch_len; ++i) {
                        tradeList.prefetchTrades1(begin[i]);
                        tradeList.prefetchTrades2(begin[i - prefetch_len]);
                    }

                    for(size_t i = 2 * prefetch_len; i < size; ++i) {
                        tradeList.prefetchTrades1(begin[i]);
                        tradeList.prefetchTrades2(begin[i - prefetch_len]);
                        update<false>(u, begin[i - 2 * prefetch_len]);
                    }

                    for(size_t i = size; i < size + prefetch_len; ++i) {
                        tradeList.prefetchTrades2(begin[i - prefetch_len]);
                        update<false>(u, begin[i - 2 * prefetch_len]);
                    }

                    for(size_t i = size + prefetch_len; i < size + 2 * prefetch_len; ++i)
                        update<false>(u, begin[i - 2 * prefetch_len]);
                };


                // Send u's neighbours
                if (edge_between_uv) update<Directed>(u, v);
                prefetched_update(u, disjoint_neighbours.data(), u_setsize);

                // Send v's neighbours
                if (edge_between_vu) update<Directed>(v, u);
                prefetched_update(v, disjoint_neighbours.data() + u_setsize, disjoint_neighbours.size() - u_setsize);

                // Send common neighbours
                for (const auto common : common_neighbours) {
                    update<Directed>(u, common);
                    update<Directed>(v, common);
                }
            }
        }
    }

    template <bool Directed>
	void loadFromGraph(const trade_vector &trades) {
        maxDegree = G.maxDegree();

        tradeList.initialize(trades);

        // Insert to adjacency list, directed according trades
        // TODO: make parallel
        G.forEdges([&](node u, node v) {update<Directed>(u, v);});
    }

    template <bool Directed>
	void restructureGraph(const trade_vector &trades) {
        // TODO: replace by forEdge
        nodepair_vector edges = adjList.getEdges();

        adjList.restructure();
        tradeList.initialize(trades);

        for (const auto edge : edges) {update<Directed>(edge.first, edge.second);}
    }

    template <bool Directed>
    void update(const node a, const node b) {
	    if (Directed) {
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

};


} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_CURVEBALL_IMPL_H
