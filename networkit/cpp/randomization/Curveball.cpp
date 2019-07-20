/*
 * Curveball.cpp
 *
 *  Created on: 26.05.2018
 *      Author:  Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <networkit/auxiliary/Random.hpp>

#include <networkit/randomization/Curveball.hpp>
#include "curveball/CurveballImpl.hpp"

namespace NetworKit {

Curveball::Curveball(const Graph &G, bool allowSelfLoops, bool isBipartite) :
    impl(new CurveballDetails::CurveballImpl{G, allowSelfLoops, isBipartite})
{
    if (allowSelfLoops && !G.isDirected()) {
        throw std::runtime_error("Self loops are only supported for directed graphs");
    }

    if (!allowSelfLoops && G.numberOfSelfLoops()) {
        throw  std::runtime_error("Self loops are forbidden but input graph contains some");
    }

    if (G.isWeighted()) {
        throw std::runtime_error("GlobalCurveball supports only unweighted graphs");
    }
}

// We have to define a "default" destructor here, since the definition of
// CurveballDetails::CurveballImpl is not known in the header file
Curveball::~Curveball() = default;

void Curveball::run(const CurveballDetails::trade_vector& trades) {
    impl->run(trades);
}

void Curveball::run(CurveballUniformTradeGenerator& generator) {
    impl->run(generator);
}

Graph Curveball::getGraph(bool) const {
    return impl->getGraph();
}

std::string Curveball::toString() const  {
    return "Curveball";
}

count Curveball::getNumberOfAffectedEdges() const {
    return impl->getNumberOfAffectedEdges();
}

} // namespace NetworKit
