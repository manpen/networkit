/*
 * DegreeSequenceAGWGenerator.hpp
 *
 *  Created on: 28.02.2014
 *      Author: cls
 */

#ifndef DEGREE_SEQUENCE_AGW_GENERATOR_H
#define DEGREE_SEQUENCE_AGW_GENERATOR_H

#include <vector>

#include <tlx/container/btree_set.hpp>

#include <networkit/graph/Graph.hpp>
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
    count num_nodes;            ///< number of nodes
    count sum_degrees;          ///< sum of degrees
    count num_edges;            ///< number of edges, i.e. m=M/2
    count max_degree;           ///< maximum entry in degrees
    count num_two_paths;        ///< number of two path sum_i d_i * (d_i - 1)
    count num_max_loops;        ///< maximal number of loops
    count num_max_double_edges; ///< maximal number of double edges
    std::vector<node> node_stubs;

    /**
     * Computes all statistics and values based on the degree sequences alone.
     * These results remain valid even after a rejection in later phases and won't be repeated.
     */
    void initializeStats();

    // representation of multigraph which get iteratively merged into simple graph
    std::vector<Edge>    edges;
    tlx::btree_set<Edge> edge_set;
    std::vector<node>    nodes_with_loop;
    std::vector<Edge>    double_edges;
    std::vector<count>   simple_degrees;

    count num_simple_twopaths; ///< Number of simple 2-path, i.e. (u,v,w) with (u,v) and (v,w) being simple edges and v without loop

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
     *
     * @return true: continue, false: reject
     */
    bool sampleInitialMultiGraph();

    bool removeLoops();
    bool removeDoubleEdges();

    /**
     * Create a graph from the vector edges assuming it satisfies the degree distribution in seq.
     */
    Graph graphFromEdgeList() const;

    static constexpr count two_path_from_degree(count deg) {
        // don't care about underflows, since if (x-1) underflows x will be zero and so will be x(x-1) ...
        // cast only to avoid warnings
        return static_cast<count>(deg * (static_cast<int64_t>(deg)-1));
    }
};

} // namespace NetworKit

#endif // DEGREE_SEQUENCE_AGW_GENERATOR_H
