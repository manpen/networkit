/*
 * HyperbolicGeneratorQuadTree.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu),
 *      Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef HYPERBOLICGENERATOR_QUADTREE_H_
#define HYPERBOLICGENERATOR_QUADTREE_H_

#include <vector>
#include <networkit/generators/HyperbolicGenerator.hpp>

namespace NetworKit {

/**
 * Random Hyperbolic Graph Generator with support for positive temperature as described in
 * "Querying Probabilistic Neighborhoods in Spatial Data Sets Efficiently" by Moritz von Looz
 * and Henning Meyerhenke, presented at IWOCA 2016.
 *
 * @warning If you are not interested in this specific algorithm but simply the fastest,
 * consider using @a HyperbolicGenerator which heuristically selects "the best" implementation.
 */
class HyperbolicGeneratorQuadTree final : public Hyperbolic::GeneratorBase {
public:
    using Hyperbolic::GeneratorBase::GeneratorBase; // pull all constructors into public scope

    /**
     * Set the capacity of a quadtree leaf.
     *
     * @param capacity Tuning parameter, default value is 1000 for T == 0 and 10 for T > 0
     */
    void setLeafCapacity(count capacity) noexcept {
        capacity = capacity;
    }

    /**
     * When using a theoretically optimal split, the quadtree will be flatter, but running time usually longer.
     * @param split Whether to use the theoretically optimal split. Defaults to false
     */
    void setTheoreticalSplit(bool split) noexcept {
        theoreticalSplit = split;
    }

    void setBalance(double balance) noexcept {
        balance = balance;
    }

private:
    /**
     * tuning parameters
     */
    count capacity{0}; //= 0 mean automatic
    bool theoreticalSplit{false};
    double balance{0.5};

    Graph generateGraph() override;
};

} // namespace NetworKit

#endif /* HYPERBOLICGENERATOR_QUADTREE_H_ */
