/*
 * DegreeSequenceAGWGenerator.hpp
 *
 *  Created on: 28.02.2014
 *      Author: cls
 */

#ifndef DEGREE_SEQUENCE_AGW_GENERATOR_H
#define DEGREE_SEQUENCE_AGW_GENERATOR_H

#include <vector>

#include <networkit/Graph.hpp>
#include <networkit/generators/StaticDegreeSequenceGenerator.hpp>

namespace NetworKit {

/**
 * Implementation of the INC-GEN / INC-POWERLAW as described in
 * "Fast uniform generation of random graphs with given degree
 * sequences" by Aramn, Gao, Wormald (arXiv: 1905.03446)
 */
class DegreeSequenceAGWGenerator : public StaticDegreeSequenceGenerator {
public:
    DegreeSequenceAGWGenerator() = delete;
    explicit DegreeSequenceAGWGenerator(const std::vector<count> &degrees);

	virtual Graph generate() override;

private:
    // statistics of and values derived from sequence
    const count num_nodes;            ///< number of nodes
    const count sum_degrees;          ///< sum of degrees
    const count num_edges;            ///< number of edges, i.e. m=M/2
    const count max_degree;           ///< maximum entry in degrees
    const count num_two_paths;        ///< number of two path sum_i d_i * (d_i - 1)
    const count num_max_loops;        ///< maximal number of loops
    const count num_max_double_edges; ///< maximal number of double edges

    std::vector<Node> node_stubs;
    std::vector<Edge> edges;

    std::vector<node> nodes_with_loop;
    std::vector<Edge> double_edges;

    /**
     * Samples the edge list of a multigraph subject to the following constraints:
     *  - no node has two or more loops
     *  - there are at most @param allowed_num_loops many nodes with a loop
     *  - there are no edges with multiplicity of three or more
     *  - there are no more than @param allowed_num_double_edges double edges
     *
     * The result is stored in the vectors edges, nodes_with_loop, double_edges.
     * The edge list edges will be sorted and edges {u, v} will be store as tuples with
     * u < v.
     *
     * Precondition: degrees must be realizable.
     */
    void sampleInitialMultiGraph(count allowed_num_loops, count allowed_num_double_edges);

    void removeLoops();

    /**
     * Create a graph from the vector edges assuming it satisfies the degree distribution in seq.
     */
    Graph graphFromEdgeList() const;
};

} // namespace NetworKit

#endif // DEGREE_SEQUENCE_AGW_GENERATOR_H
