/*
 * HyperbolicGenerator.cpp
 *
 *      Authors: Mustafa Özdayi and Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 *
 * This generator contains algorithms described in two publications.
 *
 * For T=0, the relevant publication is
 * "Generating massive complex networks with hyperbolic geometry faster in practice" by
 * Moritz von Looz, Mustafa Özdayi, Sören Laue and Henning Meyerhenke, presented at HPEC 2016.
 *
 * For T>0, it is
 * "Querying Probabilistic Neighborhoods in Spatial Data Sets Efficiently" by Moritz von Looz
 * and Henning Meyerhenke, presented at IWOCA 2016.
 *
 * The model of hyperbolic random graphs is presented in
 * "Hyperbolic geometry of complex networks. Physical Review E, 82:036106, Sep 2010." by
 *   Dmitri Krioukov, Fragkiskos Papadopoulos, Maksim Kitsak, Amin Vahdat, and Marian Boguna
 *
 */

#include <cmath>

#include <cstdlib>
#include <random>
#include <assert.h>
#include <omp.h>
#include <algorithm>

#include "../graph/GraphBuilder.h"
#include "HyperbolicGenerator.h"
#include "quadtree/Quadtree.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {

/**
 * Construct a generator for n nodes and m edges
 */
HyperbolicGenerator::HyperbolicGenerator(count n, double avgDegree, double plexp, double T) {
	nodeCount = n;
	if (plexp <= 2) throw std::runtime_error("Exponent of power-law degree distribution must be > 2");
	if (T < 0 || T == 1) throw std::runtime_error("Temperature must be non-negative and not 1.");//Really necessary? Graphs with T=1 can be generated, only their degree is not controllable
	if (avgDegree >= n) throw std::runtime_error("Average Degree must be at most n-1");
	if (T < 1) {
		alpha = 0.5*(plexp-1);
	} else {
		alpha = 0.5*(plexp-1)/T;
	}

	R = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha, T);
	std::cout << "Chose radius R " << R << std::endl << std::flush;
	temperature=T;
	initialize();
}

void HyperbolicGenerator::initialize() {
	if (temperature == 0) {
		capacity = 1000;
	} else {
		capacity = 10;
	}
	theoreticalSplit = false;
	threadtimers.resize(omp_get_max_threads());
	balance = 0.5;
}

Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, R, alpha, temperature);
}

Graph HyperbolicGenerator::generate(count n, double R, double alpha, double T) {
	assert(R > 0);
	vector<double> angles(n);
	vector<double> radii(n);

	//sample points randomly
	HyperbolicSpace::fillPoints(angles, radii, R, alpha);
	vector<index> permutation(n);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	//can probably be parallelized easily, but doesn't bring much benefit
	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> anglecopy(n);
	vector<double> radiicopy(n);

	#pragma omp parallel for
	for (omp_index j = 0; j < static_cast<omp_index>(n); j++) {
		anglecopy[j] = angles[permutation[j]];
		radiicopy[j] = radii[permutation[j]];
	}

	INFO("Generated Points");
	return generate(anglecopy, radiicopy, R, T);
}

Graph HyperbolicGenerator::generateCold(const vector<double> &angles, const vector<double> &radii, double R) {
	const count n = angles.size();
	assert(radii.size() == n);

	for (index i = 0; i < n; i++) {
		assert(radii[i] < R);
	}

	vector<index> permutation(n);
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
		permutation[i] = i;
	}

	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> bandRadii = getBandRadii(n, R);
	//Initialize empty bands
	vector<vector<Point2D<double>>> bands(bandRadii.size() - 1);
	//Put points to bands
	#pragma omp parallel for
	for (omp_index j = 0; j < static_cast<omp_index>(bands.size()); j++){
		for (index i = 0; i < n; i++){
			double alias = permutation[i];
			if (radii[alias] >= bandRadii[j] && radii[alias] <= bandRadii[j+1]){
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
	vector<vector<double>> bandAngles(bandCount);
	#pragma omp parallel for
	for (omp_index i=0; i < static_cast<omp_index>(bandCount); i++){
		const count currentBandSize = bands[i].size();
		bandAngles[i].resize(currentBandSize);
		for(index j=0; j < currentBandSize; j++) {
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
	vector<double> empty;
	GraphBuilder result(n, false, false);

	#pragma omp parallel
	{
		index id = omp_get_thread_num();
		threadtimers[id].start();
		#pragma omp for schedule(guided) nowait
		for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
			const double coshr = cosh(radii[i]);
			const double sinhr = sinh(radii[i]);
			count expectedDegree = (4/PI)*n*exp(-(radii[i])/2);
			vector<index> near;
			near.reserve(expectedDegree*1.1);
			Point2D<double> pointV(angles[i], radii[i], i);
			for(index j = 0; j < bandCount; j++){
				if(directSwap || bandRadii[j+1] > radii[i]){
					double minTheta, maxTheta;
					std::tie (minTheta, maxTheta) = getMinMaxTheta(angles[i], radii[i], bandRadii[j], R);
					//minTheta = 0;
					//maxTheta = 2*PI;
					vector<Point2D<double>> neighborCandidates = getPointsWithinAngles(minTheta, maxTheta, bands[j], bandAngles[j]);

					const count sSize = neighborCandidates.size();
					for(index w = 0; w < sSize; w++){
						double deltaPhi = PI - abs(PI-abs(angles[i] - neighborCandidates[w].getX()));
						if (coshr*cosh(neighborCandidates[w].getY())-sinhr*sinh(neighborCandidates[w].getY())*cos(deltaPhi) <= coshR) {
							if (neighborCandidates[w].getIndex() != i){
								near.push_back(neighborCandidates[w].getIndex());
							}
						}
					}
				}
			}
			if (directSwap) {
				auto newend = std::remove(near.begin(), near.end(), i); //no self loops!
				if (newend != near.end()) {
					assert(newend+1 == near.end());
					assert(*(newend)==i);
					near.pop_back();//std::remove doesn't remove element but swaps it to the end
				}
				result.swapNeighborhood(i, near, empty, false);
			} else {
				for (index j : near) {
					if (j >= n) ERROR("Node ", j, " prospective neighbour of ", i, " does not actually exist. Oops.");
					if(radii[j] > radii[i] || (radii[j] == radii[i] && angles[j] < angles[i]))
						result.addHalfEdge(i,j);
				}
			}
		}
		threadtimers[id].stop();
	}
	timer.stop();
	INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
	return result.toGraph(!directSwap, true);
}

Graph HyperbolicGenerator::generate(const vector<double> &angles, const vector<double> &radii, double R, double T) {
	if (T < 0) throw std::runtime_error("Temperature cannot be negative.");
	if (T == 0) return generateCold(angles, radii, R);
	assert(T > 0);

	Aux::Timer timer;
	timer.start();
	const count n = angles.size();
	assert(radii.size() == n);

	for (index i = 0; i < n; i++) {
		assert(radii[i] < R);
	}

	assert(alpha > 0);

	vector<index> permutation(n);
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); i++) {
		permutation[i] = i;
	}

	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> bandRadii = getBandRadii(n, R, 0.97);

	vector<double> bandLimitCosh(bandRadii.size());
	vector<double> bandLimitSinh(bandRadii.size());

	for (index i = 0; i < bandRadii.size(); i++) {
		bandLimitCosh[i] = cosh(bandRadii[i]);
		bandLimitSinh[i] = sinh(bandRadii[i]);
	}

	//Initialize empty bands
	vector<vector<Point2D<double>>> bands(bandRadii.size() - 1);
	//Put points to bands
	count pushed = 0;
	#pragma omp parallel for reduction(+:pushed)
	for (omp_index j = 0; j < static_cast<omp_index>(bands.size()); j++){
		for (index i = 0; i < n; i++){
			double alias = permutation[i];
			if (radii[alias] >= bandRadii[j] && radii[alias] < bandRadii[j+1]){
				bands[j].push_back(Point2D<double>(angles[alias], radii[alias], alias));
				pushed++;
			}
		}
	}
	if (pushed != n) throw std::runtime_error("Bands filled inconsistently: " + std::to_string(pushed) + " != " + std::to_string(n));

	Aux::Timer bandTimer;
	bandTimer.start();
	const count bandCount = bands.size();

	//1.Extract band angles to use them later, can create a band class to handle this more elegantly
	vector<vector<double>> bandAngles(bandCount);
	vector<vector<double>> bandCoshR(bandCount);
	vector<vector<double>> bandSinhR(bandCount);
	#pragma omp parallel for
	for (omp_index j=0; j < static_cast<omp_index>(bandCount); j++){
		const count currentBandSize = bands[j].size();
		bandAngles[j].resize(currentBandSize);
		bandCoshR[j].resize(currentBandSize);
		bandSinhR[j].resize(currentBandSize);
		for(index i=0; i < currentBandSize; i++) {
			bandAngles[j][i] = bands[j][i].getX();
			bandCoshR[j][i] = cosh(bands[j][i].getY());
			bandSinhR[j][i] = sinh(bands[j][i].getY());
		}
		if (!std::is_sorted(bandAngles[j].begin(), bandAngles[j].end())) {
			throw std::runtime_error("Points in bands must be sorted.");
		}
	}
	bandTimer.stop();
	INFO("Extracting band angles took ", bandTimer.elapsedMilliseconds(), " milliseconds.");

	const double coshR = cosh(R);
	assert(radii.size() == n);

	//now define lambda
	double beta = 1/T;
	assert(!std::isnan(beta));
	auto edgeProb = [beta, R](double distance) -> double {return 1 / (exp(beta*(distance-R)/2)+1);};

	auto leftOf = [](double phi, double psi){

	if (phi < PI) {
		return psi > phi && psi < phi + PI;
	} else {
		return psi > phi || psi < phi - PI;
	}

	};

	auto angleDist = [](double phi, double psi){ return PI - std::abs(PI-std::abs(phi - psi)); };

	//get Graph
	GraphBuilder result(n, false, false);//no direct swap with probabilistic graphs
	count totalCandidates = 0;

	for (index bandIndex = 0; bandIndex < bandCount; bandIndex++) {
		const omp_index bandSize = static_cast<omp_index>(bands[bandIndex].size());
		#pragma omp parallel for reduction(+:totalCandidates)
		for (omp_index bandSweepIndex = 0; bandSweepIndex < bandSize; bandSweepIndex++) {
			index i = bands[bandIndex][bandSweepIndex].getIndex();

			const double coshRI = bandCoshR[bandIndex][bandSweepIndex];
			const double sinhRI = bandSinhR[bandIndex][bandSweepIndex];

			double mirrorphi;
			if (angles[i] >= PI) mirrorphi = angles[i] - PI;
			else mirrorphi = angles[i] + PI;

			for(index j = bandIndex; j < bandCount; j++){
				const double& coshBandR = bandLimitCosh[j+1];
				const double& sinhBandR = bandLimitSinh[j+1];
				const double& coshBandRLower = bandLimitCosh[j];
				const double& sinhBandRLower = bandLimitSinh[j];

				if (bandAngles[j].size() == 0) {
					continue;
				}

				assert(bandAngles[j].size() > 0);
				//get point in b_j with next angle clockwise
				const auto it = std::lower_bound(bandAngles[j].begin(), bandAngles[j].end(), angles[i]);
				const int nextBandIndex = std::distance(bandAngles[j].begin(), it);
				int cIndex = nextBandIndex;
				assert(cIndex >= 0);

				double upperBoundProb = 1;

				auto confirmPoint = [&](int cIndex){
					//if (bands[j][cIndex].getIndex() == i) {
					//	return;
					//}
					totalCandidates += 1;

					double deltaPhi = angleDist(angles[i], bandAngles[j][cIndex]);
					double coshDist = coshRI*bandCoshR[j][cIndex]-sinhRI*bandSinhR[j][cIndex]*cos(deltaPhi);
					double distance;
					if (coshDist >= 1) distance = acosh(coshDist);
					else distance = 0;

					//double distance = HyperbolicSpace::nativeDistance(angles[i], radii[i], bandAngles[j][cIndex], );
					double q = edgeProb(distance);
					q = q / upperBoundProb; //since the candidate was selected by the jumping process, we have to adjust the probabilities
					if (q > 1) {
						double candidateR = bands[j][cIndex].getY();
						throw std::runtime_error("Upper bound " + std::to_string(upperBoundProb) + " was wrong: ("
								+ std::to_string(angles[i]) + ", " + std::to_string(radii[i]) + "), ("
								+ std::to_string(bandAngles[j][cIndex]) + ", "+ std::to_string(candidateR) + ")"
								+ ", real distance " + std::to_string(distance)
								+ ", real probability " + std::to_string(edgeProb(distance))
								+ ", deltaPhi " + std::to_string(deltaPhi));
					}
					assert(q <= 1);
					assert(q >= 0);

					//accept?
					double acc = Aux::Random::real();
					if (acc < q) {
						result.addHalfEdge(i, bands[j][cIndex].getIndex());
					}
				};

				auto advanceIndex = [&](int cIndex){
					//advance! - careful, the following is only an approximation - but should be right
					const double deltaPhi = angleDist(angles[i], bandAngles[j][cIndex]);

					const double coshDist = coshRI*coshBandR-sinhRI*sinhBandR*cos(deltaPhi);
					const double coshDistLower = coshRI*coshBandRLower-sinhRI*sinhBandRLower*cos(deltaPhi);
					const double epsilon = 0.001;//to avoid issues caused by rounding errors

					const double lowerBoundDistance = std::max(0.0, std::min(acosh(coshDistLower),  acosh(coshDist)-(bandRadii[j+1]-bandRadii[j]))  - epsilon);
					if (lowerBoundDistance <= R) {
						upperBoundProb = 1;
						return 0.0;
					}

//					const double candidateR = bands[j][cIndex].getY();
//					const double coshDistNeighbor = coshRI*bandCoshR[j][cIndex]-sinhRI*bandSinhR[j][cIndex]*cos(deltaPhi);
//					const double distNeighbor = acosh(coshDistNeighbor);
//
//					if (distNeighbor < lowerBoundDistance) {
//						throw std::runtime_error("Trigonometry error! Lower bound " + std::to_string(lowerBoundDistance)
//						+ ", derived from " + std::to_string(acosh(coshDist)) + " - " + std::to_string((bandRadii[j+1]-bandRadii[j]))
//						+ ", but neighbor distance is " + std::to_string(distNeighbor)
//						+ ". Neighbor at (" + std::to_string(bandAngles[j][cIndex]) + ", " + std::to_string(candidateR) + ")"
//						+ ", bands at radii " + std::to_string(bandRadii[j]) + " and " + std::to_string(bandRadii[j+1])
//						+ ". Query point at (" + std::to_string(angles[i]) + ", " + std::to_string(radii[i]) + "). "
//						+ "Distance to higher band " + std::to_string(HyperbolicSpace::nativeDistance(angles[i], radii[i], bandAngles[j][cIndex], bandRadii[j+1])) + ". "
//						+ "Distance to lower band " + std::to_string(HyperbolicSpace::nativeDistance(angles[i], radii[i], bandAngles[j][cIndex], bandRadii[j])) + ". "
//						+ "Distance between bands " + std::to_string(HyperbolicSpace::nativeDistance(bandAngles[j][cIndex], bandRadii[j], bandAngles[j][cIndex], bandRadii[j+1])) + ". "
//						+ "Native distance to neighbor " + std::to_string(HyperbolicSpace::nativeDistance(angles[i], radii[i], bandAngles[j][cIndex], candidateR))
//						);
//					}

					upperBoundProb = std::min(edgeProb(lowerBoundDistance)*1.01, 1.0);

					double probdenom = std::log(1-upperBoundProb);
					double random = Aux::Random::real();
					double delta = std::log(random) / probdenom;
					return delta;
				};

				int pointsSkipped = 0;
				while (cIndex < bandAngles[j].size() && pointsSkipped < bandAngles[j].size() && leftOf(angles[i], bandAngles[j][cIndex])) {
					//add point or not
					if (bands[j][cIndex].getY() >= radii[i]) {
						confirmPoint(cIndex);
					}

					double delta = advanceIndex(cIndex);

					cIndex += int(delta) + 1;
					pointsSkipped += int(delta) + 1;

					if (cIndex >= bandAngles[j].size()) {
						cIndex -= bandAngles[j].size();
					}
				}

				cIndex = nextBandIndex - 1;
				if (cIndex < 0) {
					cIndex += bandAngles[j].size();
				}

				upperBoundProb = 1;
				pointsSkipped = 0;
				while (cIndex >= 0 && cIndex < bandAngles[j].size() && pointsSkipped < bandAngles[j].size() && leftOf(bandAngles[j][cIndex], angles[i]) ) {
					//add point or not
					if (bands[j][cIndex].getY() >= radii[i]) {
						confirmPoint(cIndex);
					}

					double delta = advanceIndex(cIndex);
					pointsSkipped += int(delta) + 1;

					cIndex -= int(delta) + 1;
					if (cIndex < 0) {
						cIndex += bandAngles[j].size();
					}
				}
			}
		}
	}
	std::cout << "Candidates tested: " << totalCandidates << std::endl << std::flush;
	return result.toGraph(true, true);

}
}
