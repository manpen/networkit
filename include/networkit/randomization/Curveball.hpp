/*
 * Curveball.h
 *
 *  Created on: 26.05.2018
 *      Author:  Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */
#ifndef RANDOMIZATION_CURVEBALL_H
#define RANDOMIZATION_CURVEBALL_H

#include <memory>
#include <utility>

#include <networkit/Globals.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

// forward declaration for pImpl
namespace CurveballDetails { class CurveballImpl; }

class Curveball : public Algorithm {
public:

	explicit Curveball(const Graph &G, bool allowSelfLoops = false, bool isBipartite = false);

	virtual ~Curveball();

	void run() override final {
		throw std::runtime_error("run() is not supported by this algorithm; use run(trades)");
	};

	void run(const std::vector<std::pair<node, node> >& trades);

	///! Returns a copy of the randomized graph.
	///! The @a parallel flag is ignored and will eventually be removed
	Graph getGraph(bool parallel = true) const;

	std::string toString() const final;

	bool isParallel() const final {
		return false;
	}

	count getNumberOfAffectedEdges() const;

private:
	std::unique_ptr<CurveballDetails::CurveballImpl> impl;
};

} // ! namespace NetworKit

#endif // ! RANDOMIZATION_CURVEBALL_H
