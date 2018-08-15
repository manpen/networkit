#include <gtest/gtest.h>

#include "../FastGraphBuilder.h"
#include "../Graph.h"

#include <set>
#include <map>
#include <random>
#include <utility>

namespace NetworKit {

template <typename FGB>
class FastGraphBuilderGTest: public testing::Test {};

using BuilderTypesToTest = ::testing::Types<FastGraphBuilder<false>, FastGraphBuilder<true> >;
TYPED_TEST_CASE(FastGraphBuilderGTest, BuilderTypesToTest);

TYPED_TEST(FastGraphBuilderGTest, UndirectedClique) {
    using builder_type = TypeParam;

    const node n = 100;
    builder_type builder(n, builder_type::supportWeights, false);

    size_t num_edges = 0;
    #pragma omp parallel for schedule(dynamic, 100), reduction(+:num_edges)
    for(node i=1; i<n; i++) {
        const auto tid = omp_get_thread_num();
        for(node j=0; j<i; j++) {
            builder.addEdge(tid, i, j, i * j);
            num_edges++;
        }
    }

    Graph G = builder.toGraph();

    ASSERT_EQ(G.numberOfNodes(), n);
    ASSERT_EQ(G.numberOfEdges(), n * (n-1) / 2);

    size_t count = 0;
    G.forEdges([&count] (node u, node v, edgeweight w) {
        count++;
        if (builder_type::supportWeights) {
            ASSERT_DOUBLE_EQ(w, u * v);
        }
    });

    ASSERT_EQ(count, G.numberOfEdges());
}

TYPED_TEST(FastGraphBuilderGTest, UndirectedErdosReyni) {
    std::default_random_engine prng{1};
    for(size_t iter=0; iter<1; iter++) {
        using builder_type = TypeParam;
        using edge = std::pair<node, node>;
        const node n = 10000;
        const edgeid m = 2 * n;

        bool directed = false;


// generate m random edges in a random permutation
        std::vector<std::pair<edge, edgeweight> > edges_to_insert;
        {
            std::map<edge, edgeweight> map;

            std::uniform_int_distribution<node> ndistr(0, n - 1);
            std::uniform_real_distribution<edgeweight> wdistr;

            while (map.size() < m) {
                const node u = ndistr(prng);
                const node v = ndistr(prng);
                const edgeweight w = builder_type::supportWeights ? wdistr(prng) : 1.0;

                const edge e = directed || (u > v) ? edge{u, v} : edge{v, u};

                map[e] = w;
            }

            edges_to_insert.reserve(map.size());
            for (const auto &p : map) {
                edges_to_insert.emplace_back(p.first, p.second);
            }

            std::shuffle(edges_to_insert.begin(), edges_to_insert.end(), prng);
        }

// feed to graph builder
        const int threads = 129;
        builder_type builder(n, builder_type::supportWeights, directed, threads);

        {
            std::uniform_int_distribution<int> tdistr(0, threads - 1);
            for (const auto &p : edges_to_insert)
                builder.addEdge(tdistr(prng), p.first.first, p.first.second, p.second);
        }

// build graph
        Graph G = builder.toGraph();
        ASSERT_EQ(G.numberOfNodes(), n);
        ASSERT_EQ(G.numberOfEdges(), edges_to_insert.size());

        std::vector<std::pair<edge, edgeweight> > edges_in_graph;
        edges_in_graph.reserve((1 + directed) * edges_to_insert.size());

        std::sort(edges_to_insert.begin(), edges_to_insert.end());

        G.forEdges([&](node u, node v, edgeweight w) {
            edges_in_graph.emplace_back(edge{u, v}, w);

            edge e{u,v};
            auto it = std::lower_bound(edges_to_insert.cbegin(), edges_to_insert.cend(), e,
                                       [] (const std::pair<edge, edgeweight>& i, const edge& e) {return i.first < e;});
            ASSERT_NE(it, edges_to_insert.cend());
            ASSERT_EQ(it->first, e);
        });

// compare
        std::sort(edges_in_graph.begin(), edges_in_graph.end());

        ASSERT_GE(edges_in_graph.size(), edges_to_insert.size());

        auto git = edges_in_graph.cbegin();
        size_t i = 0;
        for (auto iit = edges_to_insert.cbegin(); iit != edges_to_insert.cend(); ++iit, ++git, ++i) {

            ASSERT_EQ(iit->first, git->first) << i;
            ASSERT_DOUBLE_EQ(iit->second, git->second) << i;

            // selfloops are emitted twice by NetworKit, but only contained once in the input stream
            if (git->first.first == git->first.second) {
                ASSERT_NE(git + 1, edges_in_graph.cend());
                ASSERT_EQ(git->first, (git + 1)->first);
                ++git;
            }
        }
    }
}




} // namespace NetworKit
