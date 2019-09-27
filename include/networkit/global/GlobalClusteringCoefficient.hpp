/*
 * GlobalClusteringCoefficient.h
 *
 *  Created on: 12.11.2013
 */

#ifndef GLOBALCLUSTERINGCOEFFICIENT_H_
#define GLOBALCLUSTERINGCOEFFICIENT_H_

#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup global
 */
class GlobalClusteringCoefficient {

public:  
	double approximate(const Graph& G, count k);
};

} /* namespace NetworKit */
#endif /* GLOBALCLUSTERINGCOEFFICIENT_H_ */
