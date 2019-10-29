/*
 * DGSWriter.h
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#ifndef DGSWRITER_H_
#define DGSWRITER_H_

#include <string>
#include <vector>

#include <networkit/dynamics/GraphEvent.hpp>

namespace NetworKit {

/**
 * @ingroup dynamics
 */
class DGSWriter {
public:
    DGSWriter() = default;

    void write(std::vector<GraphEvent>& stream, const std::string& path);

};

} /* namespace NetworKit */

#endif /* DGSWRITER_H_ */
