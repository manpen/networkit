/*
 * DissimilarityMeasure.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "DissimilarityMeasure.h"

namespace NetworKit {


double DissimilarityMeasure::getDissimilarity(const Graph&, const Cover&, const Cover&) {
	throw std::runtime_error("Not implemented");
}

} /* namespace NetworKit */
