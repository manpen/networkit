/*
 * DynamicPubWebGenerator.h
 *
 *  Created on: 15.01.2014
 *      Author: Henning
 */
// networkit-format

#ifndef DYNAMICPUBWEBGENERATOR_H_
#define DYNAMICPUBWEBGENERATOR_H_

#include <map>
#include <vector>

#include <networkit/dynamics/GraphEvent.hpp>
#include <networkit/generators/DynamicGraphGenerator.hpp>
#include <networkit/generators/PubWebGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class DynamicPubWebGenerator : public DynamicGraphGenerator {
public:
    DynamicPubWebGenerator(count numNodes, count numberOfDenseAreas, coordinate neighborhoodRadius,
                           count maxNumberOfNeighbors, bool writeInitialGraphToStream = true);

    Graph getGraph() const { return G; }

    /**
     * Generate event stream.
     *
     * @param[in]	nSteps	number of time steps in the event stream
     */
    std::vector<GraphEvent> generate(count nSteps) override;

    /// Returns a map of coordinates that were updated.
    const std::map<node, Point2D> &getNewCoordinates() const { return newCoordinates; }
    std::map<node, Point2D> moveNewCoordinates() { return std::move(newCoordinates); }

    /// Returns a vector of the currently valid coordinates
    const std::vector<Point2D> &getCoordinates() const { return coordinates; }
    // no moveCoordinates, as generator needs its own copy for the next run of generate!

protected:
    PubWebGenerator initGen; // multiple inheritance did not work with different generate functions
    std::map<node, Point2D> newCoordinates; //<! new and changed coordinates
    std::vector<Point2D> coordinates;       //<! vector of all coordinates
    bool writeInitialGraphToStream;         // if true, on first call, write initial graph to stream
};

} /* namespace NetworKit */
#endif /* DYNAMICPUBWEBGENERATOR_H_ */
