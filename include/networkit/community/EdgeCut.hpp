/*
 * EdgeCut.h
 *
 *  Created on: Jun 20, 2013
 *      Author: Henning
 */

#ifndef EDGECUT_H_
#define EDGECUT_H_

#include <networkit/community/QualityMeasure.hpp>

namespace NetworKit {

/**
 * @ingroup community
 */
class EdgeCut: public QualityMeasure {
public:
    virtual double getQuality(const Partition& zeta, const Graph& G);
};

} /* namespace NetworKit */
#endif /* EDGECUT_H_ */
