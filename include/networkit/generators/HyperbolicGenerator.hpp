/*
 * HyperbolicGenerator.h
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu), Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef HYPERBOLICGENERATOR_H_
#define HYPERBOLICGENERATOR_H_

#include <vector>

#include <networkit/graph/Graph.hpp>
#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * Generator for Random Hyperbolic Graphs
 *
 * This class is a dispatch which tries to identify the most suited implementation
 * for the requested set of parameters.
 */
class HyperbolicGenerator: public StaticGraphGenerator {
public:
	/**
	 * @param[in] n Number of nodes
	 * @param[in] k Target average degree
	 * @param[in] plexp Target exponent of power-law distribution
	 * @param[in] T Temperature
	 */
	HyperbolicGenerator(count n=10000, double avgDegree=6, double plexp=3, double T=0);

	/**
	 * @param[in] angles Angular coordinates of points (have to be in the interval [0:2Pi))
	 * @param[in] radii  Polar coordinates of points (have to be in the interval [0:R))
	 * @param[in] R Target radius of hyperbolic plane
	 * @param[in] plexp Target exponent of power-law distribution (here used only as a tuning parameter)
	 * @param[in] T Temperature
	 */
	HyperbolicGenerator(std::vector<double> angles, std::vector<double> radii, double R, double plexp, double T);

	/// @return the expected powerlaw exponent based on the alpha parameter and the temperature T
	static constexpr double alphaToPLE(double alpha, double T) {
		return (T < 1) ? 2.0 * alpha + 1   : 2.0 * T * alpha + 1;
	}

	/// @return the alpha parameter based on the requested powerlaw exponent plexp and the temperature T
	static constexpr double PLEToAlpha(double plexp, double T) {
		return (T < 1) ? 0.5 * (plexp - 1) : 0.5 * (plexp - 1) / T;
	}

	/// @return Graph to be generated according to parameters specified in constructor freeing memory by deleting the input point set.
	Graph generate() override {
		return generate(false);
	}

	/// @return Graph to be generated according to parameters specified in constructor keeping the input point set.
	Graph generateKeepingInput() {
		return generate(true);
	}

	/// @return Model Parameter Alpha (i.e., the displacement parameter)
	double getAlpha() const noexcept {
		return alpha;
	}

	/// @return Model Parameter T (i.e., the temperature)
	double getT() const noexcept {
		return temperature;
	}

	/// @return Model Paramter R (i.e., the target radius)
	double getR() const noexcept {
		return R;
	}

	/**
	 * @return  Angles used to generate the graph
	 * @warning The data is destroyed if generate is called with keep_input = false (default).
	 */
	const std::vector<double>& getAngles() const noexcept {
		return angles;
	}

	std::vector<double>& getAngles() noexcept {
		return angles;
	}

	/**
	 * @return  Weights used to generate the graph
	 * @warning The data is destroyed if generate is called with keep_input = false (default).
	 */
	const std::vector<double>& getRadii() const noexcept {
		return radii;
	}

	std::vector<double>& getRadii() noexcept {
		return radii;
	}

	/// @return Expected powerlaw exponent of the generated graph
	double getPowerlawExponent() const noexcept {
		return plexp;
	}

private:
	Graph generate(bool keepInput);

	/**
	 * graph parameters
	 */
	count nodeCount;
	double R;
	double alpha;
	double plexp;
	double temperature;

	std::vector<double> angles;
	std::vector<double> radii;
};

} // namespace NetworKit

#endif /* HYPERBOLICGENERATOR_H_ */
