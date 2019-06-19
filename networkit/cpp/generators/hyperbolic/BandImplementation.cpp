/*
 * HyperbolicGenerator.cpp
 *
 *      Authors: Mustafa zdayi, Moritz v. Looz (moritz.looz-corswarem@kit.edu), Manuel Penschuck (networkit@manuel.jetzt)
 */

#include <cmath>

#include <cstdlib>
#include <random>
#include <assert.h>
#include <omp.h>
#include <algorithm>

#include <networkit/graph/GraphBuilder.hpp>
#include <networkit/generators/hyperbolic/BandImplementation.hpp>
#include <networkit/generators/quadtree/Quadtree.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/Parallel.hpp>

namespace NetworKit {
namespace HyperbolicGenerators {

Graph BandImplementation::generate(const std::vector<double> &angles, const std::vector<double> &radii, double R) {
	const count n = angles.size();
	assert(radii.size() == n);

	for (index i = 0; i < n; i++) {
		assert(radii[i] < R);
	}

	threadtimers.resize(omp_get_max_threads());

	std::vector<index> permutation(n);
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
		permutation[i] = i;
	}

	//can probably be parallelized easily, but doesn't bring much benefit
	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles, &radii](index i, index j) {
		return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);
	});

	std::vector<double> bandRadii = getBandRadii(n, R);
	//Initialize empty bands
	std::vector<std::vector<Point2D<double>>> bands(bandRadii.size() - 1);
	//Put points to bands
	#pragma omp parallel for
	for (omp_index j = 0; j < static_cast<omp_index>(bands.size()); j++) {
		for (index i = 0; i < n; i++) {
			double alias = permutation[i];
			if (radii[alias] >= bandRadii[j] && radii[alias] <= bandRadii[j + 1]) {
				bands[j].push_back(Point2D<double>(angles[alias], radii[alias], alias));
			}
		}
	}

	const count bandCount = bands.size();
	const double coshR = cosh(R);
	assert(radii.size() == n);

	Aux::Timer bandTimer;
	bandTimer.start();

	//1.Extract band angles to use them later, can create a band class to handle this more elegantly
	std::vector<std::vector<double>> bandAngles(bandCount);
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(bandCount); i++) {
		const count currentBandSize = bands[i].size();
		bandAngles[i].resize(currentBandSize);
		for (index j = 0; j < currentBandSize; j++) {
			bandAngles[i][j] = bands[i][j].getX();
		}
		if (!std::is_sorted(bandAngles[i].begin(), bandAngles[i].end())) {
			throw std::runtime_error("Points in bands must be sorted.");
		}
	}
	bandTimer.stop();
	INFO("Extracting band angles took ", bandTimer.elapsedMilliseconds(), " milliseconds.");

	//2.Insert edges
	Aux::Timer timer;
	timer.start();
	std::vector<double> empty;
	GraphBuilder result(n, false, false);

	#pragma omp parallel
	{
		index id = omp_get_thread_num();
		threadtimers[id].start();
		#pragma omp for schedule(guided) nowait
		for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
			const double coshr = cosh(radii[i]);
			const double sinhr = sinh(radii[i]);
			count expectedDegree = (4 / PI) * n * exp(-(radii[i]) / 2);
			std::vector<index> near;
			near.reserve(expectedDegree * 1.1);
			Point2D<double> pointV(angles[i], radii[i], i);
			for (index j = 0; j < bandCount; j++) {
				if (directSwap || bandRadii[j + 1] > radii[i]) {
					double minTheta, maxTheta;
					std::tie(minTheta, maxTheta) = getMinMaxTheta(angles[i], radii[i], bandRadii[j], R);
					//minTheta = 0;
					//maxTheta = 2*PI;
					std::vector<Point2D<double>> neighborCandidates = getPointsWithinAngles(minTheta, maxTheta, bands[j], bandAngles[j]);

					const count sSize = neighborCandidates.size();
					for (index w = 0; w < sSize; w++) {
						double deltaPhi = PI - std::abs(PI - std::abs(angles[i] - neighborCandidates[w].getX()));
						if (coshr * std::cosh(neighborCandidates[w].getY()) - sinhr * sinh(neighborCandidates[w].getY()) * std::cos(deltaPhi) <=
							coshR) {
							if (neighborCandidates[w].getIndex() != i) {
								near.push_back(neighborCandidates[w].getIndex());
							}
						}
					}
				}
			}
			if (directSwap) {
				auto newend = std::remove(near.begin(), near.end(), i); //no self loops!
				if (newend != near.end()) {
					assert(newend + 1 == near.end());
					assert(*(newend) == i);
					near.pop_back();//std::remove doesn't remove element but swaps it to the end
				}
				result.swapNeighborhood(i, near, empty, false);
			} else {
				for (index j : near) {
					if (j >= n) ERROR("Node ", j, " prospective neighbour of ", i, " does not actually exist. Oops.");
					if (radii[j] > radii[i] || (radii[j] == radii[i] && angles[j] < angles[i]))
						result.addHalfEdge(i, j);
				}
			}
		}
		threadtimers[id].stop();
	}
	timer.stop();
	INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
	return result.toGraph(!directSwap, true);
}

} // namespace HyperbolicGenerators
} // namespace NetworKit
