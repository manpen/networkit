/*
 * GraphBuilderBenchmark.cpp
 *
 *  Created on: 04.12.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com), Manuel Penschuck <networkit@manuel.jetzt>
 *
 */

#include <benchmark/benchmark.h>
#include "../../generators/ErdosRenyiEnumerator.h" // this is a header-only dependency

#include "../GraphBuilder.h"
#include "../FastGraphBuilder.h"


namespace NetworKit {

class GraphBuilderBenchmark : public ::benchmark::Fixture {
public:
	static void param_scan(benchmark::internal::Benchmark* b) {
		for (int nodes = 1llu << 10; nodes <= 1llu << 18; nodes *= 16)
			for (int degree = 1; degree <= 256; degree *= 16)
				for (int directed = 0; directed <= 1; directed++)
					b->Args({nodes, degree, directed});
	}

protected:
	ErdosRenyiEnumerator<> getEnumerator(const benchmark::State &state) {
		const node nNodes = state.range(0);
		const edgeid nEdges = nNodes * state.range(1);
		const bool directed = state.range(2);
		const double prob = 1.0 * nEdges / (directed ? nNodes * nNodes : nNodes * (nNodes - 1) / 2);
		return {nNodes, prob, directed};
	}
};


BENCHMARK_DEFINE_F(GraphBuilderBenchmark, EnumeratorBaseline)(benchmark::State &state) {
	auto ere = this->getEnumerator(state);

	count num_edges = 0;
	for (auto _ : state) {
		num_edges += ere.forEdgesParallel([&](int tid, node u, node v) {
			benchmark::DoNotOptimize(u);
			benchmark::DoNotOptimize(v);
		});
	}

	state.SetItemsProcessed(num_edges);
}

BENCHMARK_REGISTER_F(GraphBuilderBenchmark, EnumeratorBaseline)->Apply(GraphBuilderBenchmark::param_scan);

BENCHMARK_DEFINE_F(GraphBuilderBenchmark, GraphBuilderParFillOnly)(benchmark::State &state) {
    auto ere = this->getEnumerator(state);

    count num_edges = 0;
    for (auto _ : state) {
        GraphBuilder builder(state.range(0));

        num_edges += ere.forEdgesParallel([&](int tid, node u, node v) {
            builder.addHalfEdge(u, v);
        });

        benchmark::DoNotOptimize(builder); // should be unnecessary, but let's be on the safe side
    }

    state.SetItemsProcessed(num_edges);
}
BENCHMARK_REGISTER_F(GraphBuilderBenchmark, GraphBuilderParFillOnly)->Apply(GraphBuilderBenchmark::param_scan);

BENCHMARK_DEFINE_F(GraphBuilderBenchmark, GraphBuilderSeqBuild)(benchmark::State &state) {
    auto ere = this->getEnumerator(state);

    count num_edges = 0;
    for (auto _ : state) {
        GraphBuilder builder(state.range(0));

        num_edges += ere.forEdgesParallel([&](int tid, node u, node v) {
            builder.addHalfEdge(u, v);
        });

        auto G = builder.toGraph(true, false);
    }

    state.SetItemsProcessed(num_edges);
}
BENCHMARK_REGISTER_F(GraphBuilderBenchmark, GraphBuilderSeqBuild)->Apply(GraphBuilderBenchmark::param_scan);

BENCHMARK_DEFINE_F(GraphBuilderBenchmark, GraphBuilderParBuild)(benchmark::State &state) {
    auto ere = this->getEnumerator(state);

    count num_edges = 0;
    for (auto _ : state) {
        GraphBuilder builder(state.range(0));

        num_edges += ere.forEdgesParallel([&](int tid, node u, node v) {
            builder.addHalfEdge(u, v);
        });

        auto G = builder.toGraph(true, true);
    }

    state.SetItemsProcessed(num_edges);
}
BENCHMARK_REGISTER_F(GraphBuilderBenchmark, GraphBuilderParBuild)->Apply(GraphBuilderBenchmark::param_scan);

BENCHMARK_DEFINE_F(GraphBuilderBenchmark, GraphBuilderFastBuilderUnweighted)(benchmark::State &state) {
    auto ere = this->getEnumerator(state);

    count num_edges = 0;
    for (auto _ : state) {
        // test unweighted gbp
        FastGraphBuilder<false> builder(state.range(0), false, state.range(2));

        num_edges += ere.forEdgesParallel([&](int tid, node u, node v) {
            builder.addEdge(tid, u, v);
        });

        auto G = builder.toGraph();

        std::cout << "\n\n\n";
    }

    state.SetItemsProcessed(num_edges);
}
BENCHMARK_REGISTER_F(GraphBuilderBenchmark, GraphBuilderFastBuilderUnweighted)->Apply(GraphBuilderBenchmark::param_scan);

BENCHMARK_DEFINE_F(GraphBuilderBenchmark, GraphBuilderFastBuilderWeighted)(benchmark::State &state) {
    auto ere = this->getEnumerator(state);

    count num_edges = 0;
    for (auto _ : state) {
        // test unweighted gbp
        FastGraphBuilder<false> builder(state.range(0), true, state.range(2));

        num_edges += ere.forEdgesParallel([&](int tid, node u, node v) {
            builder.addEdge(tid, u, v, 0.5);
        });

        auto G = builder.toGraph();

        std::cout << "\n\n\n";
    }

    state.SetItemsProcessed(num_edges);
}
BENCHMARK_REGISTER_F(GraphBuilderBenchmark, GraphBuilderFastBuilderWeighted)->Apply(GraphBuilderBenchmark::param_scan);


} // ! namespace NetworKit