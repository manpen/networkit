/*
 * HyperbolicGeneratorQuadTree.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu), Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef HYPERBOLICGENERATOR_QUADTREE_H_
#define HYPERBOLICGENERATOR_QUADTREE_H_

#include <vector>

namespace NetworKit {
namespace HyperbolicGenerators {
/**
 * Random Hyperbolic Graph Generator with support for positive temperature as described in
 * "Querying Probabilistic Neighborhoods in Spatial Data Sets Efficiently" by Moritz von Looz
 * and Henning Meyerhenke, presented at IWOCA 2016.
 *
 * @warning This is an implementation class only which does not satisfy the StaticGenerator interface
 * by design. Only use if you know what you are doing; prefer usage of Dispatcher class HyperbolicGenerator!
 */
class QuadTreeImplementation {
public:
	QuadTreeImplementation() = default;

	/**
	 * @param[in] angles Pointer to angles of node positions
	 * @param[in] radii Pointer to radii of node positions
	 * @param[in] r radius of poincare disk to place nodes in
	 * @param[in] thresholdDistance Edges are added for nodes closer to each other than this threshold
	 * @return Graph to be generated according to parameters
	 */
	Graph generate(const std::vector<double> &angles, const std::vector<double> &radii, double R, double T, double alpha);

	/**
	 * Set the capacity of a quadtree leaf.
	 *
	 * @param capacity Tuning parameter, default value is 1000 for T == 0 and 10 for T > 0
	 */
	void setLeafCapacity(count capacity) {
		capacity = capacity;
	}

	/**
	 * When using a theoretically optimal split, the quadtree will be flatter, but running time usually longer.
	 * @param split Whether to use the theoretically optimal split. Defaults to false
	 */
	void setTheoreticalSplit(bool split) {
		theoreticalSplit = split;
	}

	void setBalance(double balance) {
		balance = balance;
	}

private:
	/**
	 * tuning parameters
	 */
	count capacity{0}; //= 0 mean automatic
	bool theoreticalSplit{false};
	double balance{0.5};
};

} // namespace HyperbolicGenerators
} // namespace NetworKit

#endif /* HYPERBOLICGENERATOR_QUADTREE_H_ */
