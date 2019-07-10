/*
 * EdgeListDistance.hpp
 *
 *  Created on: 06.07.2019
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef EDGE_LIST_DIFFERENCE_HPP_
#define EDGE_LIST_DIFFERENCE_HPP_

#include <limits>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Distance measures between edges lists.
 *
 * Given two graphs G_1 = (V_1, E_2) and G_2 = (V_2, E_2), this algorithm
 * computes the number of common and disjoint edges, as well as
 * the jaccard index, the overlapping coefficient, and the perturbation score
 * on the edge sets.
 *
 * The algorithm has a runtime of O(m + dmax log dmax) and a space complexity
 * of O(dmax * P) where m is the number of edges, P the number of threads and
 * dmax the maximal degree.
 */
class EdgeSetDistance : public Algorithm {
public:
    EdgeSetDistance(const Graph &g1, const Graph &g2);

// Algorithm overloads
    void run() override;

    std::string toString() const override {
        return "EdgeListDistance";
    }

    bool isParallel() const override {
        return true;
    }

// Scores
    /**
     * Size of the intersection of both edge sets, i.e. |intersection(E_1, E_2)|.
     * @warning Requires previous call of run()
     */
    count sizeOfIntersection() const {
        assureFinished();
        return num_common_edges;
    }

    /**
     * Size of the union of both edge sets, i.e. |union(E_1, E_2)|.
     * @warning Requires previous call of run()
     */
    count sizeOfUnion() const {
        assureFinished();
        return g1.numberOfEdges() + g2.numberOfEdges() - num_common_edges;
    }

    /**
     * The jaccard index is given as |intersection(E_1, E_2)| / |E_1 union E_2|.
     * @warning Requires previous call of run()
     */
    double jaccardIndex() const {
        assureFinished();
        return static_cast<double>(sizeOfIntersection()) / sizeOfUnion();
    }

    /**
     * The overlap coefficient (also known as Szymkiewicz–Simpson coefficient)
     * is given as |intersection(E_1, E_2)| / min(|E_1|, |E_2|).
     * @warning Requires previous call of run().
     */
    double overlapCoefficient() const {
        assureFinished();
        return static_cast<double>(sizeOfIntersection()) / std::min(g1.numberOfEdges(), g2.numberOfEdges());
    }

    /**
     * The perturbation score is given as |intersection(E_1, E_2)| / |E_1|.
     * @warning The perturbation score is only defined if |E_1| = |E_2| and is then
     * identical to the overlap_coefficient.
     * @warning Requires previous call of run().
     */
    double perturbationScore() const {
        assureFinished();
        if (g1.numberOfEdges() != g2.numberOfEdges())
            throw std::runtime_error("Perturbation score only defined if the number of edges match in both graphs.");

        return overlapCoefficient();
    }

private:
    const Graph &g1;
    const Graph &g2;

    count num_common_edges;  ///< size of the intersection of E_1 and E_2

};

} // namespace NetworKit

#endif // EDGE_SET_DIFFERENCE_HPP_
