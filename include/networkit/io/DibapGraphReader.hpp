/*
 * DibapGraphReader.h
 *
 *  Created on: Jun 12, 2013
 *      Author: Henning
 */
// networkit-format

#ifndef DIBAPGRAPHREADER_H_
#define DIBAPGRAPHREADER_H_

#ifndef NETWORKIT_WINDOWS

#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>
#include <networkit/viz/Point.hpp>

namespace NetworKit {

/**
 * @ingroup io
 * TODO: class documentation
 */
class DibapGraphReader final : public GraphReader {
public:
    DibapGraphReader() = default;

    Graph read(const std::string &path) override;

    const std::vector<Point<coordinate>> &getCoordinates() const { return coordinates; }

    std::vector<Point<coordinate>> moveCoordinates() { return std::move(coordinates); }

private:
    std::vector<Point<coordinate>> coordinates;
};

} /* namespace NetworKit */

#endif // NETWORKIT_WINDOWS

#endif /* DIBAPGRAPHREADER_H_ */
