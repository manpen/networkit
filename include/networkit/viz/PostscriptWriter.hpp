/*
 * PostscriptWriter.h
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#ifndef POSTSCRIPTWRITER_H_
#define POSTSCRIPTWRITER_H_

#include <string>
#include <fstream>
#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/viz/Point.hpp>

namespace NetworKit {

/**
 * @ingroup viz
 * EPS output of graphs with 2D coordinates
 */
class PostscriptWriter {
public:
    /**
     * @param[in] isTorus Specifies whether the visualization square is treated as torus,
     * i.e. with wrap-around boundaries (edge can leave the square and enter at the opposite
     * side. By default, it is set to false.
     */
    PostscriptWriter(bool isTorus = false);

    /**
     * Outputs an EPS file with name @a filename of the graph @a g with 2D coordinates.
     * The colors are chosen to visualize the specified @a clustering.
     * @param[in] g Graph to be visualized.
     * @param[in] clustering Clustering of the graph, visualized by different colors.
     * @param[in] filename Name of file to write to.
     */
    void write(const Graph& g, const std::vector<coord2d>& coordinates, const Partition& clustering, const std::string& filename);

    /**
     * Outputs an EPS file with name @a filename of the graph @a g with 2D coordinates.
     * @param[in] g Graph to be visualized.
     * @param[in] filename Name of file to write to.
     */
    void write(const Graph& g, const std::vector<coord2d>& coordinates, const std::string& filename);

private:
    bool wrapAround;

    coord2d ps_size;
    coord2d ps_border;
    coord2d ps_min;
    coord2d ps_max;

    void init(std::string filename, std::ofstream& file);
    void writeHeader(std::ofstream& file);
    void writeMacros(std::ofstream& file);
    void writeClustering(const Graph& g, const std::vector<coord2d>& coordinates, const Partition& clustering, std::ofstream& file);

};

} /* namespace NetworKit */
#endif /* POSTSCRIPTWRITER_H_ */
