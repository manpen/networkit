/*
 * HyperbolicGeneratorQuadTree.cpp
 *
 *      Authors: Mustafa zdayi, Moritz v. Looz (moritz.looz-corswarem@kit.edu), Manuel Penschuck (networkit@manuel.jetzt)
 *
 */

#include <cmath>

#include <cstdlib>
#include <random>
#include <assert.h>
#include <omp.h>
#include <algorithm>

#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/generators/HyperbolicGeneratorQuadTree.hpp>
#include <networkit/generators/quadtree/Quadtree.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/auxiliary/Parallel.hpp>

namespace NetworKit {

Graph HyperbolicGeneratorQuadTree::generateGraph() {
    if (temperature < 0)
        throw std::runtime_error("Temperature cannot be negative.");

    if (capacity == 0)
        capacity = temperature > 0 ? 10 : 1000;

    /**
     * fill Quadtree
     */
    Aux::Timer timer;
    timer.start();
    index n = angles.size();
    assert(radii.size() == n);

    assert(alpha > 0);
    Quadtree<index, false> quad(R, theoreticalSplit, alpha, capacity, balance);

    for (index i = 0; i < n; i++) {
        assert(0 <= radii[i] && radii[i] < R);
        assert(0 <= angles[i] && angles[i] < 2*PI);
        quad.addContent(i, angles[i], radii[i]);
    }

    quad.trim();
    timer.stop();
    INFO("Filled Quadtree, took ", timer.elapsedMilliseconds(), " milliseconds.");

    assert(quad.size() == n);

    const bool anglesSorted = std::is_sorted(angles.begin(), angles.end());

    //now define lambda
    const double beta = 1.0 / temperature;
    assert(beta == beta);
    auto edgeProb = [beta, this](double distance) -> double { return 1 / (exp(beta * (distance - R) / 2) + 1); };

    //get Graph
    GraphBuilder result(n, false, false);//no direct swap with probabilistic graphs
    count totalCandidates = 0;
    #pragma omp parallel for
    for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
        std::vector<index> near;
        totalCandidates += quad.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), edgeProb, anglesSorted, near);
        for (index j : near) {
            if (j >= n) ERROR("Node ", j, " prospective neighbour of ", i, " does not actually exist. Oops.");
            if (j > i) {
                result.addHalfEdge(i, j);
            }
        }

    }

    DEBUG("Candidates tested: ", totalCandidates);
    return result.toGraph(true, true);
}

} // namespace NetworKit
