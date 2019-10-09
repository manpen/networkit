/*
 * PartitionReader.h
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef PARTITIONREADER_H_
#define PARTITIONREADER_H_

#include <fstream>

#include <networkit/structures/Partition.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * @ingroup io
 */
class PartitionReader {

public:
    /**
     * Read a clustering from a file. File format:
     * 		line n contains cluster id of node (n - 1)
     *
     * @param[in]	path	Path to file.
     */
    virtual Partition read(std::string path);
};

} /* namespace NetworKit */
#endif /* PARTITIONREADER_H_ */
