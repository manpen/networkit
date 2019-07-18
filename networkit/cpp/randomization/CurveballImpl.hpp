/*
 * CurveballImpl.h
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
#ifndef RANDOMIZATION_CURVEBALL_IMPL_H
#define RANDOMIZATION_CURVEBALL_IMPL_H

#include <cassert>
#include <utility>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

#include "curveball/AdjacencyList.hpp"

namespace NetworKit {
namespace CurveballDetails {

// Global Definitions
using trade_descriptor = std::pair<node, node>;
using tradeid = node;
using count = node;

using trade_vector = std::vector<trade_descriptor>;
using nodepair_vector = std::vector<std::pair<node, node>>;

constexpr node INVALID_NODE = std::numeric_limits<node>::max();
constexpr count LISTROW_END = std::numeric_limits<count>::max();
constexpr tradeid TRADELIST_END = std::numeric_limits<tradeid>::max();


class CurveballMaterialization {

protected:
	const CurveballAdjacencyList &adjacencyList;

public:
	CurveballMaterialization(const CurveballAdjacencyList &adj_list);

	Graph toGraph(bool parallel);

protected:
	void toGraphParallel(Graph &G);
	void toGraphSequential(Graph &G);
};

class TradeList {
public:
	using edge_vector = std::vector<std::pair<node, node>>;
	using offset_vector = std::vector<tradeid>;
	using tradeid_vector = std::vector<tradeid>;
	using trade = trade_descriptor;
	using trade_vector = std::vector<trade>;
	using tradeid_it = tradeid_vector::const_iterator;

protected:
	tradeid_vector tradeList;
	offset_vector offsets;
	const node numNodes;

public:
	TradeList(const node num_nodes);

	// Receives the edge_vector to initialize
	TradeList(const trade_vector &trades, const node num_nodes);

	// Initialize method
	void initialize(const trade_vector &trades);

	// No Copy Constructor
	TradeList(const TradeList &) = delete;

	tradeid_it getTrades(const node nodeid) const {
		assert(nodeid < numNodes);

		return tradeList.begin() + offsets[nodeid];
	}

	void incrementOffset(const node nodeid) {
		assert(nodeid < numNodes);
		assert(1 <= offsets[nodeid + 1] - offsets[nodeid]);

		offsets[nodeid]++;
	}

	node numberOfNodes() const { return numNodes; }
};

class CurveballIM {
public:
	CurveballIM(const Graph &G);

	void run(const trade_vector &trades);

	count getNumberOfAffectedEdges() const {
		assert(hasRun);
		return numAffectedEdges;
	}

	Graph getGraph(bool parallel) const;

	nodepair_vector getEdges() const;

protected:
	const Graph &G;
	const node numNodes;

	bool hasRun;
	AdjacencyList adjList;
	TradeList tradeList;
	count maxDegree;
	edgeid numAffectedEdges; // affected half-edges

	void loadFromGraph(const trade_vector &trades);

	void restructureGraph(const trade_vector &trades);

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
