/*
 * EdgeScoreBlender.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann
 */

#ifndef EDGESCOREBLENDER_H
#define EDGESCOREBLENDER_H

#include <networkit/graph/Graph.hpp>
#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

class EdgeScoreBlender : public EdgeScore<double> {

public:

    EdgeScoreBlender(const Graph &G, const std::vector<double> &attribute0, const std::vector<double> &attribute1, const std::vector<bool> &selection);

    virtual double score(edgeid eid) override;
    virtual double score(node u, node v) override;
    virtual void run() override;

private:
    const std::vector<double> &attribute0, &attribute1;
    const std::vector<bool> &selection;
};

} // namespace NetworKit

#endif // EDGESCOREBLENDER_H
