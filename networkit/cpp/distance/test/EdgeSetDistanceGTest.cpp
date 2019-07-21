#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <random>

#include <networkit/distance/EdgeSetDistance.hpp>

namespace NetworKit {

class EdgeSetDistanceGTest : public testing::Test {
protected:
    // This test performs random edits (adding of edges, deleting of nodes
    // and its incident edges, restoring isolated nodes) independently to two graphs.
    // The changes are also tracked in edge sets
    template <bool Directed>
    void runTest() {
        const node n = 30;
        const int numRounds = 8*n*n;
        using Edge = std::pair<node, node>;

        std::array<Graph, 2> graphs{ Graph(n, false, Directed), Graph(n, false, Directed) };
        std::array<std::set<Edge>, 2> edges;

        std::mt19937_64 gen;

        std::uniform_int_distribution<int> act_distr(0, 1);
        std::uniform_int_distribution<node> node_distr(0, n - 1);

        auto normalise = [] (node u, node v)  {
            return Directed ? Edge{u,v} : static_cast<Edge>(std::minmax(u, v));
        };

        for(int i = 0; i < numRounds; ++i) {
            const auto act = std::uniform_int_distribution<int>(0, 1)(gen);
            auto& act_graph = graphs[act];

            const auto action = std::uniform_real_distribution<double>{}(gen);

            if (action < 0.5) {
                // Add edge
                Edge e;
                int retries = 10;

                bool validEdge = true;

                do {
                    e = normalise(node_distr(gen), node_distr(gen));

                    if (!retries--) {
                        validEdge = false;
                        break;
                    }
                } while (
                    !act_graph.hasNode(e.first) ||
                    !act_graph.hasNode(e.second) ||
                    !edges[act].insert(e).second); // non existent edge between two existing node

                if(validEdge)
                    act_graph.addEdge(e.first, e.second);

            } else if (action < 0.9) {
                if (!act_graph.numberOfEdges())
                    continue;

                auto e =  *edges[act].cbegin();
                act_graph.removeEdge(e.first, e.second);
                edges[act].erase(edges[act].begin());

            } else {
                // Delete or restore node
                const auto v = node_distr(gen);

                if (act_graph.hasNode(v)) {
                    // Remove all incident edges from our set copy
                    act_graph.forNeighborsOf(v, [&] (node vv) {
                        edges[act].erase(normalise(v, vv));
                    });

                    if (Directed) {
                        act_graph.forInNeighborsOf(v, [&] (node vv) {
                            edges[act].erase({vv, v});
                        });
                    }

                    act_graph.removeNode(v);

                } else {
                    act_graph.restoreNode(v);
                }


            }

            ASSERT_EQ(graphs[0].numberOfEdges(), edges[0].size());
            ASSERT_EQ(graphs[1].numberOfEdges(), edges[1].size());

            if (i % 32)
                continue;

            count size_intersection = 0;
            for(const Edge e : edges[0]) {
                const bool in_both = (edges[1].find(e) != edges[1].end());
                size_intersection += in_both;
            }

            for(int j = 0; j < 2; ++j) { // Check for symmetry
                EdgeSetDistance alg(graphs[j], graphs[1-j]);
                alg.run();
                ASSERT_EQ(alg.sizeOfIntersection(), size_intersection);
                ASSERT_EQ(alg.sizeOfUnion(), edges[0].size() + edges[1].size() - size_intersection);
            }
        }

    }


};

TEST_F(EdgeSetDistanceGTest, testUndirectedRandom) {
    this->runTest<false>();
}

TEST_F(EdgeSetDistanceGTest, testDirectedRandom) {
    this->runTest<true>();
}

}
