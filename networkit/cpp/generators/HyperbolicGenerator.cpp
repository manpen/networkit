/*
 * HyperbolicGenerator.cpp
 *
 *      Authors: Manuel Penschuck (networkit@manuel.jetzt)
 *
 */

#include <cmath>

#include <assert.h>
#include <algorithm>

#include <networkit/auxiliary/Parallel.hpp>

#include <networkit/geometric/HyperbolicSpace.hpp>

#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/generators/hyperbolic/BandImplementation.hpp>
#include <networkit/generators/hyperbolic/QuadTreeImplementation.hpp>


namespace NetworKit {

/**
 * Construct a generator for n nodes and m edges
 */
HyperbolicGenerator::HyperbolicGenerator(count n, double avgDegree, double plexp, double T) :
	nodeCount(n),
	alpha(PLEToAlpha(plexp, T)),
	plexp(plexp),
	temperature(T)
{
	if (plexp <= 2)
		throw std::runtime_error("Exponent of power-law degree distribution must be > 2");

	if (T < 0 || T == 1)
		throw std::runtime_error("Temperature must be non-negative and not 1.");

	if (avgDegree >= n)
		throw std::runtime_error("Average Degree must be at most n-1");

	R = HyperbolicSpace::getTargetRadius(n, n * avgDegree / 2, alpha, T);
}

HyperbolicGenerator::HyperbolicGenerator(std::vector<double> angles, std::vector<double> radii, double R, double plexp, double T) :
	nodeCount(angles.size()),
	R(R),
	alpha(PLEToAlpha(plexp, T)),
	plexp(plexp),
	temperature(T),
	angles(std::move(angles)),
	radii(std::move(radii))
{
	if (plexp <= 2)
		throw std::runtime_error("Exponent of power-law degree distribution must be > 2");

	if (T < 0 || T == 1)
		throw std::runtime_error("Temperature must be non-negative and not 1.");

	if (this->angles.size() != this->radii.size() || !this->angles.size())
		throw std::runtime_error("Number of angles and radii has to match and to be positive");
}

Graph HyperbolicGenerator::generate(bool keepInput) {
	assert(R > 0);
	if (angles.empty()) {
		angles.resize(nodeCount);
		radii.resize(nodeCount);

		// need to generate points
		HyperbolicSpace::fillPoints(angles, radii, R, alpha);
	}

	INFO("Generated Points");

	auto clear_input = [&] {
		if (keepInput)
			return;

		angles.clear(); angles.shrink_to_fit();
		radii.clear(); radii.shrink_to_fit();
	};

	if (temperature == 0.0) {
		auto G = HyperbolicGenerators::BandImplementation().generate(angles, radii, R);
		clear_input();
		return G;

	} else {
		//sample points randomly
		std::vector<index> permutation(nodeCount);

		index p = 0;
		std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

		//can probably be parallelized easily, but doesn't bring much benefit
		Aux::Parallel::sort(permutation.begin(), permutation.end(), [&] (index i, index j) {
			return std::tie(angles[i], radii[i]) < std::tie(angles[j], radii[j]);});

		std::vector<double> anglecopy(nodeCount);
		std::vector<double> radiicopy(nodeCount);

		#pragma omp parallel for
		for (omp_index j = 0; j < static_cast<omp_index>(nodeCount); j++) {
			anglecopy[j] = angles[permutation[j]];
			radiicopy[j] = radii[permutation[j]];
		}

		clear_input();
		return HyperbolicGenerators::QuadTreeImplementation().generate(anglecopy, radiicopy, R, temperature, alpha);
	}
}

} // namespace NetworKit
