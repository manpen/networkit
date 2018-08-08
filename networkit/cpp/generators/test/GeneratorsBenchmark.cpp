/*
 * GeneratorsBenchmark.cpp
 *
 *  Created on: May 29, 2013
 *      Author: forigem
 */

#ifndef NOGTEST

#include <omp.h>
#include <random>
#include <functional>

#include "GeneratorsBenchmark.h"
#include "../../auxiliary/Log.h"
#include "../../auxiliary/Parallel.h"
#include "../../auxiliary/Parallelism.h"

#include "../ErdosRenyiEnumerator.h"
#include "../HyperbolicGenerator.h"
#include "../DynamicHyperbolicGenerator.h"
#include "../BarabasiAlbertGenerator.h"
#include "../ChungLuGenerator.h"
#include "../MocnikGenerator.h"
#include "../MocnikGeneratorBasic.h"
#include "../../graph/GraphBuilder.h"

namespace NetworKit {

constexpr node graphBuilderNodes = 100000;
constexpr double graphBuilderProb = 0.001;
constexpr double graphBuilderExpEdges = 0.5 * graphBuilderNodes * (graphBuilderNodes + 1) * graphBuilderProb;

TEST_F(GeneratorsBenchmark, benchmarkGraphBuilderParFillSeqBuild) {
	GraphBuilder builder(graphBuilderNodes);
	count m_actual = 0;
	auto t1 = timeOnce([&]() {
		ErdosRenyiEnumerator ere(graphBuilderNodes, graphBuilderProb, false);
		ere.forEdgesParallel([&](int tid, node u, node v) {
			builder.addHalfEdge(u, v);
		});
	});
	auto t2 = timeOnce([&]() {
		auto G = builder.toGraph(true, false);
		m_actual = G.numberOfEdges();
	});
	EXPECT_NEAR(m_actual / graphBuilderExpEdges, 1.0, 0.1);
	std::cout << "parallelForNodePairs + toGraphSequentiel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";
}

TEST_F(GeneratorsBenchmark, benchmarkGraphBuilderParFillParBuild) {
	GraphBuilder builder(graphBuilderNodes);
	count m_actual = 0;
	auto t1 = timeOnce([&]() {
		ErdosRenyiEnumerator ere(graphBuilderNodes, graphBuilderProb, false);
		ere.forEdgesParallel([&](int tid, node u, node v) {
			builder.addHalfEdge(u, v);
		});
	});
	auto t2 = timeOnce([&]() {
		auto G = builder.toGraph(true, true);
		m_actual = G.numberOfEdges();
	});
	EXPECT_NEAR(m_actual / graphBuilderExpEdges, 1.0, 0.1);
	std::cout << "parallelForNodePairs + toGraphParallel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";
}


TEST_F(GeneratorsBenchmark, benchmarkGraphBuilder) {
	count m_actual;
	auto t1 = timeOnce([&]() {
		auto G = Graph(graphBuilderNodes);
		ErdosRenyiEnumerator ere(graphBuilderNodes, graphBuilderProb, false);
		ere.forEdges([&](int tid, node u, node v) {
			G.addEdge(u, v);
		});
		m_actual = G.numberOfEdges();
	});
	EXPECT_NEAR(m_actual / (double) graphBuilderExpEdges, 1.0, 0.1);
	std::cout << "forNodePairs + Graph.addEdge:\t\t\t\t" << t1 << " ms\n";
}

TEST_F(GeneratorsBenchmark, benchmarkBarabasiAlbertGenerator) {
	count k = 2;
	count nMax = 100000;
	count n0 = 2;

	BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, false);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());
}

TEST_F(GeneratorsBenchmark, benchBarabasiAlbertGeneratorBatagelj) {
	for (index i = 0; i < 10; ++i) {
		Aux::Random::setSeed(i, false);
		count n = Aux::Random::integer(100, 10000);
		count k = n / Aux::Random::integer(5, 20);
		BarabasiAlbertGenerator gen(k, n, 0);
		auto G = gen.generate();
		//EXPECT_TRUE(G.checkConsistency());
		//INFO(G.toString());
	}
}

TEST_F(GeneratorsBenchmark, benchBarabasiAlbertGenerator2) {
	for (index i = 0; i < 10; ++i) {
		Aux::Random::setSeed(i, false);
		count n = Aux::Random::integer(100, 10000);
		count k = n / Aux::Random::integer(5, 20);
		BarabasiAlbertGenerator gen(k, n, k, false);
		auto G = gen.generate();
		//EXPECT_TRUE(G.checkConsistency());
		//INFO(G.toString());
	}
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGenerator) {
	count n = 100000;
	HyperbolicGenerator gen(n);
	Graph G = gen.generate();
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGeneratorWithSortedNodes) {
	count n = 100000;
	double s = 1.0;
	double alpha = 1.0;
	double t = 1.0;
	vector<double> angles(n);
	vector<double> radii(n);
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	//sample points randomly

	HyperbolicSpace::fillPoints(angles, radii, s, alpha);
	vector<index> permutation(n);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	//can probably be parallelized easily, but doesn't bring much benefit
	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> anglecopy(n);
	vector<double> radiicopy(n);

	#pragma omp parallel for
	for (index j = 0; j < n; j++) {
		anglecopy[j] = angles[permutation[j]];
		radiicopy[j] = radii[permutation[j]];
	}

	Graph G = HyperbolicGenerator().generate(anglecopy, radiicopy, r, R*t);
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkDynamicHyperbolicGeneratorOnNodeMovement) {
	const count runs = 100;
	const double fractionStep = 0.01;
	const count stepCount = 100;
	const count n = 1000000;
	const double k = 6;
	const double exp = 3;
	const double moveDistance = 0.1;

	for (index i = 0; i < stepCount; i++) {
		double moveFraction = fractionStep * i;
		for (index j = 0; j < runs; j++) {
			DynamicHyperbolicGenerator dyngen(n, k, exp, moveFraction, moveDistance);
			dyngen.generate(1);
		}
	}
}

TEST_F(GeneratorsBenchmark, benchmarkParallelQuadtreeConstruction) {
	count n = 33554432;
	Quadtree<index> quad(n,1.0);
	EXPECT_EQ(quad.size(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkSequentialQuadtreeConstruction) {
	count n = 33554432;
	count capacity = 1000;
	double s =1;
	double alpha = 1;
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, s, alpha);

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),false,alpha,capacity);

	for (index i = 0; i < n; i++) {
		quad.addContent(i, angles[i], radii[i]);
	}
	EXPECT_EQ(quad.size(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGeneratorMechanicGraphs) {
	count n = 1000000;
	double k = 6;
	count m = n*k/2;
	HyperbolicGenerator gen(n, k, 3, 0.14);
	gen.setLeafCapacity(10);
	Graph G = gen.generate();
	EXPECT_NEAR(G.numberOfEdges(), m, m/10);
}

TEST_F(GeneratorsBenchmark, benchmarkChungLuGenerator) {
	count n = 100000;
    int maxDegree = 100;
	std::vector<count> vec;
	/* Creates a random weight list */
	for (index i = 0; i < n; i++){
	int grad = Aux::Random::integer(1, maxDegree);
		vec.push_back(grad);
	}
	ChungLuGenerator generator(vec);
	Graph G = generator.generate();
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkMocnikGenerator) {
	count dim = 3;
	count n = 1000000;
	double k = 2.6;
	
	MocnikGenerator Mocnik(dim, n, k);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());
	G = Mocnik.generate();
	EXPECT_FALSE(G.isEmpty());
	EXPECT_EQ(G.numberOfNodes(), n);
	EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 2000000);
}

TEST_F(GeneratorsBenchmark, benchmarkMocnikGeneratorBasic) {
	count dim = 3;
	count n = 10000;
	double k = 2.6;
	
	MocnikGenerator Mocnik(dim, n, k);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());
	G = Mocnik.generate();
	EXPECT_FALSE(G.isEmpty());
	EXPECT_EQ(G.numberOfNodes(), n);
	EXPECT_NEAR(G.numberOfEdges() * 1. / G.numberOfNodes(), std::pow(k, dim), 20000);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
