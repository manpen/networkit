/*
 * PostscriptWriterGTest.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include <gtest/gtest.h>

#include <vector>

#include <networkit/viz/PostscriptWriter.hpp>
#include <networkit/viz/MaxentStress.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/generators/ClusteredRandomGraphGenerator.hpp>
#include <networkit/io/PartitionWriter.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISGraphWriter.hpp>
#include <networkit/io/DibapGraphReader.hpp>
#include <networkit/generators/PubWebGenerator.hpp>
#include <networkit/auxiliary/Random.hpp>


namespace NetworKit {

class VizGTest : public testing::Test {};


TEST_F(VizGTest, testPostscriptWriterOnRandomGraph) {
    // create graph
    count n = 60;
    count numClusters = 3;
    double pin = 0.35;
    double pout = 0.05;

    ClusteredRandomGraphGenerator graphGen(n, numClusters, pin, pout);
    Graph G = graphGen.generate();
    std::vector<coord2d> coordinates(G.upperNodeIdBound());

    // create coordinates
    G.forNodes([&](node u) {
        coordinates[u] = {Aux::Random::probability(), Aux::Random::probability()};
    });


    // write graph to file
    std::string path = "output/testGraph.eps";
    PostscriptWriter psWriter;
    psWriter.write(G, coordinates, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

#ifndef NETWORKIT_WINDOWS
TEST_F(VizGTest, testPostscriptWriterOnRealGraph) {
    // read graph and coordinates from binary file
    DibapGraphReader reader;
    Graph G = reader.read("input/airfoil1.gi");
    const auto coordinates = reader.moveCoordinates();


    // write graph to file
    std::string path = "output/airfoil1.eps";
    PostscriptWriter psWriter;
    psWriter.write(G, Point<>::pointVectorToCoord2d(coordinates), path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}
#endif


static float edgeDistanceSum(Graph& G) {
    float dist = 0.0f;

    G.forEdges([&](node u, node v) {
        Point<float> p = G.getCoordinate(u) - G.getCoordinate(v);
        dist += p.length();
    });

    return dist;
}

TEST_F(VizGTest, testFRLayouter) {
    // create graph
    count n = 80;
    count numClusters = 3;
    double pin = 0.175;
    double pout = 0.005;

    ClusteredRandomGraphGenerator graphGen(n, numClusters, pin, pout);
    Graph G = graphGen.generate();
    G.initCoordinates();
    DEBUG("Number of edges: ", G.numberOfEdges());

    // draw (independent of clustering) and write again
    Point<float> bl(0.0, 0.0);
    Point<float> tr(1.0, 1.0);

    PostscriptWriter psWriter(true);
    psWriter.write(G, "output/testForceGraph.eps");

    // test edge distances
    float dist = edgeDistanceSum(G);
    float avg = dist / (float) G.numberOfEdges();
    DEBUG("avg edge length: ", avg);
    EXPECT_LE(avg, 0.25);
}

TEST_F(VizGTest, debugGraphDrawing) {
    // create graph
    METISGraphReader reader;
    Graph G = reader.read("input/lesmis.graph");

    // draw (independent of clustering) and write again
    Point<float> bl(0.0, 0.0);
    Point<float> tr(1.0, 1.0);

    PostscriptWriter psWriter2(true);
    psWriter2.write(G, "output/testLesmisFR.eps");

    // test edge distances
    float dist = edgeDistanceSum(G);
    float avg = dist / (float) G.numberOfEdges();
    INFO("avg edge length: ", avg);
    EXPECT_LE(avg, 0.25);

    PostscriptWriter psWriter4(true);
    psWriter4.write(G, "output/testLesmisMl.eps");

    // test edge distances
    dist = edgeDistanceSum(G);
    avg = dist / (float) G.numberOfEdges();
    INFO("avg edge length: ", avg);
    EXPECT_LE(avg, 0.25);
}


} /* namespace NetworKit */
