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

#include "curveball/AdjacencyList.hpp"
#include "curveball/TradeList.hpp"


namespace NetworKit {
namespace CurveballDetails {

using nodepair_vector = std::vector<std::pair<node, node>>;

constexpr node INVALID_NODE = std::numeric_limits<node>::max();
constexpr count LISTROW_END = std::numeric_limits<count>::max();


class CurveballIM {
    using neighbour_vector = std::vector<node>;

public:
	CurveballIM(const Graph &G) :
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

	void run(const trade_vector &trades) {
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

        for (const auto &trade : trades) {
            handler.assureRunning();

            // Trade partners u and v
            const node u = trade.first;
            const node v = trade.second;

            numAffectedEdges += adjList.degreeAt(u);
            numAffectedEdges += adjList.degreeAt(v);

            // Shift the tradeList offset for these two, currently was set to
            // trade_count
            tradeList.incrementOffset(u);
            tradeList.incrementOffset(v);

            // Retrieve respective neighbours
            // we return whether u has v in his neighbors or vice-versa
            auto organize_neighbors = [&](node node_x, node node_y) {
                auto pos = std::find(adjList.begin(node_x), adjList.end(node_x), node_y);
                if (pos == adjList.cend(node_x)) {
                    // element not found, sort anyway
                    std::sort(adjList.begin(node_x), adjList.end(node_x));

                    return false;
                } else {
                    // overwrite node_y's position with END
                    *pos = LISTROW_END;

                    // sort, such that node_y's position is at end - 1
                    std::sort(adjList.begin(node_x), adjList.end(node_x));

                    // overwrite with node_y again
                    *(adjList.end(node_x) - 1) = node_y;

                    return true;
                }
            };

            const bool u_share = organize_neighbors(u, v);
            const bool v_share = organize_neighbors(v, u);
            auto u_end = (u_share ? adjList.cend(u) - 1 : adjList.cend(u));
            auto v_end = (v_share ? adjList.cend(v) - 1 : adjList.cend(v));

            const bool shared = u_share || v_share;

            // both can't have each other, only inserted in one
            assert((!u_share && !v_share) || (u_share != v_share));

            // No need to keep track of direct positions
            // Get common and disjoint neighbors
            // Here sort and parallel scan
            common_neighbours.clear();
            disjoint_neighbours.clear();
            auto u_nit = adjList.cbegin(u);
            auto v_nit = adjList.cbegin(v);
            while ((u_nit != u_end) && (v_nit != v_end)) {
                assert(*u_nit != v);
                assert(*v_nit != u);
                if (*u_nit > *v_nit) {
                    disjoint_neighbours.push_back(*v_nit);
                    v_nit++;
                    continue;
                }
                if (*u_nit < *v_nit) {
                    disjoint_neighbours.push_back(*u_nit);
                    u_nit++;
                    continue;
                }
                // *u_nit == *v_nit
                {
                    common_neighbours.push_back(*u_nit);
                    u_nit++;
                    v_nit++;
                }
            }
            if (u_nit == u_end)
                disjoint_neighbours.insert(disjoint_neighbours.end(), v_nit, v_end);
            else
                disjoint_neighbours.insert(disjoint_neighbours.end(), u_nit, u_end);

            const count u_setsize = static_cast<count>(u_end - adjList.cbegin(u) -
                                                       common_neighbours.size());
            const count v_setsize = static_cast<count>(v_end - adjList.cbegin(v) -
                                                       common_neighbours.size());
            // v_setsize not necessarily needed

            // Reset fst/snd row
            adjList.resetRow(u);
            adjList.resetRow(v);

            tlx::random_bipartition_shuffle(disjoint_neighbours.begin(),
                                            disjoint_neighbours.end(), u_setsize, urng);

            // Assign first u_setsize to u and last v_setsize to v
            // if not existent then max value, and below compare goes in favor of
            // partner, if partner has no more neighbours as well then their values are
            // equal (max and equal) and tiebreaking is applied
            for (count counter = 0; counter < u_setsize; counter++) {
                const node swapped = disjoint_neighbours[counter];
                update(u, swapped);
            }
            for (count counter = u_setsize; counter < u_setsize + v_setsize;
                 counter++) {
                const node swapped = disjoint_neighbours[counter];
                update(v, swapped);
            }
            // Distribute common edges
            for (const auto common : common_neighbours) {
                update(u, common);
                update(v, common);
            }
            // Do not forget edge between u and v
            if (shared)
                update(u, v);

            trade_count++;
        }

        hasRun = true;

        return;
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
        G.forEdges([&](node u, node v) { update(u, v); });
    }

	void restructureGraph(const trade_vector &trades) {
        nodepair_vector edges = adjList.getEdges();

        adjList.restructure();
        tradeList.initialize(trades);

        for (const auto edge : edges) {
            update(edge.first, edge.second);
        }
    }

	inline void update(const node a, const node b) {
		const tradeid ta = *(tradeList.getTrades(a));
		const tradeid tb = *(tradeList.getTrades(b));
		if (ta < tb) {
			adjList.insertNeighbour(a, b);
			return;
		}

		if (ta > tb) {
			adjList.insertNeighbour(b, a);
			return;
		}
		// ta == tb
		{ adjList.insertNeighbour(a, b); }
	}
};


} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_CURVEBALL_IMPL_H
