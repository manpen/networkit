/*
 * GlobalCurveballGTest.cpp
 *
 *  Created on: 24.05.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#include <gtest/gtest.h>

#include <networkit/graph/Graph.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>
#include <networkit/generators/HyperbolicGenerator.hpp>

#include <networkit/randomization/Curveball.hpp>
#include <networkit/randomization/GlobalCurveball.hpp>
#include <networkit/randomization/CurveballUniformTradeGenerator.hpp>

namespace NetworKit {


/////////////////////// ALGOS
struct GlobalCurveballGTestGlobalCurveball {
    static Graph run(const Graph &G, bool allowSelfLoops) {
        GlobalCurveball algo(G, 5, allowSelfLoops);
        algo.run();
        return algo.getGraph();
    }
};

struct GlobalCurveballGTestCurveballUniform {
    static Graph run(const Graph &G, bool allowSelfLoops) {
        CurveballUniformTradeGenerator trades(G.numberOfNodes(), G.numberOfNodes());
        Curveball algo(G, allowSelfLoops);
        algo.run(trades.generate());
        algo.run(trades.generate());
        return algo.getGraph();
    }
};


using CurveballAlgos = ::testing::Types<GlobalCurveballGTestGlobalCurveball, GlobalCurveballGTestCurveballUniform>;
TYPED_TEST_CASE(GlobalCurveballGTest, CurveballAlgos);

template <typename Algo>
class GlobalCurveballGTest : public ::testing::Test  {
protected:

    void checkWithUndirectedGraph(Graph& G) {
        ASSERT_FALSE(G.isDirected());

        node numNodes = G.numberOfNodes();
        const count numTrades = 5;

        std::vector<node> degrees(numNodes + 1);

        // Add edge to node 0, if isolated node
        // If 0 itself is isolated, add new node and connect 0 to it
        G.forNodes([&](node u) {
            if (G.degree(u) > 0)
                degrees[u] = G.degree(u);
            else {
                if (u == 0) {
                    numNodes++;
                    G.addEdge(0, numNodes - 1);
                    degrees[0]++;
                    degrees[numNodes - 1] = 1;
                } else {
                    G.addEdge(u, 0);
                    degrees[0]++;
                    degrees[u] = 1;
                }
            }
        });

        Graph outG = Algo::run(G, false);

        // check degrees
        outG.forNodes([&](node u){
            ASSERT_EQ(degrees[u], outG.degree(u));
        });
    }

    void checkWithDirectedGraph(Graph& G, count selfLoops) {
        ASSERT_TRUE(G.isDirected());

        node numNodes = G.numberOfNodes();
        const count numTrades = 5;

        std::vector<node> degreesIn(numNodes + 1);
        std::vector<node> degreesOut(numNodes + 1);

        // Add edge to node 0, if isolated node
        // If 0 itself is isolated, add new node and connect 0 to it
        G.forNodes([&](node u) {
            if (G.degreeIn(u) > 0 || G.degreeOut(u) > 0) {
                degreesIn[u] = G.degreeIn(u);
                degreesOut[u] = G.degreeOut(u);
            } else {
                if (u == 0) {
                    numNodes++;
                    G.addEdge(0, numNodes - 1);
                    degreesOut[0]++;
                    degreesIn[numNodes - 1] = 1;
                } else {
                    G.addEdge(u, 0);
                    degreesIn[0]++;
                    degreesOut[u] = 1;
                }
            }
        });

        Graph outG = Algo::run(G, selfLoops);

        // check degrees
        outG.forNodes([&](node u){
            ASSERT_EQ(degreesIn[u], outG.degreeIn(u));
            ASSERT_EQ(degreesOut[u], outG.degreeOut(u));
        });

        if (selfLoops) {
            ASSERT_GT(outG.numberOfSelfLoops(), selfLoops);
        } else {
            ASSERT_EQ(outG.numberOfSelfLoops(), 0);
        }
    }
};

TYPED_TEST(GlobalCurveballGTest, testCBUndirectedErdosRenyi) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    ErdosRenyiGenerator generator(numNodes, 0.01);
    Graph G = generator.generate();

    this->checkWithUndirectedGraph(G);
}

TYPED_TEST(GlobalCurveballGTest, testCurveballDirectedErdosRenyi) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    ErdosRenyiGenerator generator(numNodes, 0.01, true, false);
    Graph G = generator.generate();

    this->checkWithDirectedGraph(G, 0);
}

TYPED_TEST(GlobalCurveballGTest, testCurveballDirectedErdosRenyiSelfLoops) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    double p = 0.01;
    ErdosRenyiGenerator generator(numNodes, 0.01, true, true);
    Graph G = generator.generate();

    this->checkWithDirectedGraph(G, numNodes * p / 2);
}

TYPED_TEST(GlobalCurveballGTest, testCurveballHyperbolic) {
    Aux::Random::setSeed(1, false);

    node numNodes = 1000;
    HyperbolicGenerator generator(numNodes);
    Graph G = generator.generate();

    this->checkWithUndirectedGraph(G);
}


} // namespace NetworKit
