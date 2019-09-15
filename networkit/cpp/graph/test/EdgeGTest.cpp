/*
 * EdgeGTest.cpp
 *
 *  Created on: 15.09.2019
 *      Author: Manuel Penschuck
 */

#include <string>
#include <sstream>

#include <gtest/gtest.h>

#include <networkit/graph/Edge.hpp>

namespace NetworKit {

class EdgeGTest : public testing::Test {
};

TEST_F(EdgeGTest, testEdgeBasics) {
    // construct
    ASSERT_EQ(Edge(2, 1, false).first,  2);
    ASSERT_EQ(Edge(2, 1, false).second, 1);
    ASSERT_EQ(Edge(3, 4, true).first,   3);
    ASSERT_EQ(Edge(3, 4, true).second,  4);

    // construct via pair
    ASSERT_EQ(Edge(std::pair<node, node>(2, 1), false).first, 2);
    ASSERT_EQ(Edge(std::pair<node, node>(2, 1), false).second, 1);
    ASSERT_EQ(Edge(std::pair<node, node>(2, 1), true ).first, 1);
    ASSERT_EQ(Edge(std::pair<node, node>(2, 1), true ).second, 2);

    // construct via tuple
    ASSERT_EQ(Edge(std::tuple<node, node>(2, 1), false).first, 2);
    ASSERT_EQ(Edge(std::tuple<node, node>(2, 1), false).second, 1);
    ASSERT_EQ(Edge(std::tuple<node, node>(2, 1), true ).first, 1);
    ASSERT_EQ(Edge(std::tuple<node, node>(2, 1), true ).second, 2);

    // compares
    ASSERT_EQ(Edge(1, 2), Edge(2, 1, true));
    ASSERT_NE(Edge(2, 1), Edge(2, 1, true));

    // copy and move
    {
        Edge e1(3, 4);
        Edge e2 = e1;
        Edge e3 = std::move(e1);
        ASSERT_EQ(e2, e3);
    }

    // flip
    {
        Edge e1{1,9};
        auto e2 = e1.flipped();

        ASSERT_NE(e1, e2);
        e1.flip();
        ASSERT_EQ(e1, e2);
    }

    // sort
    {
        Edge e1{6,5};
        auto e2 = e1.sorted();

        ASSERT_NE(e1, e2);
        e1.sort();
        ASSERT_EQ(e1, e2);
        e2.sort();
        ASSERT_EQ(e1, e2);
        e1 = e1.sorted();
        ASSERT_EQ(e1, e2);
    }

    // incident
    {
        Edge e1(1, 2), e2(2, 3), e3(3, 4);

        auto check = [&] {
            ASSERT_TRUE (e1.isIncident(1));
            ASSERT_TRUE (e1.isIncident(2));
            ASSERT_FALSE(e1.isIncident(3));

            ASSERT_FALSE(e2.isIncident(1));
            ASSERT_TRUE (e2.isIncident(3));
            ASSERT_TRUE (e2.isIncident(2));

            ASSERT_TRUE (e1.isIncident(e2));
            ASSERT_FALSE(e1.isIncident(e3));
        };

        check();
        e1.flip();
        check();
        e2.flip();
        check();
        e3.flip();
        check();
    }

    // loop
    {
        ASSERT_TRUE (Edge(1,1).isLoop());
        ASSERT_FALSE(Edge(1,2).isLoop());
    }

    // deconstruction
    {
        Edge e(5, 6);
        node u = 0, v = 0;
        std::tie(u, v) = e;

        ASSERT_EQ(u, 5);
        ASSERT_EQ(v, 6);
    }

    // print
    {
        std::stringstream ss;
        Edge e(5, 6);
        ss << e;

        std::string s = ss.str();

        ASSERT_NE(s.find('5'), std::string::npos);
        ASSERT_NE(s.find('6'), std::string::npos);
        ASSERT_LE(s.find('5'), s.find('6'));
    }

    // undirected equality
    {
        Edge e1(3, 6), e2(6, 3);

        ASSERT_FALSE(e1.isEqualDirected(e2));
        ASSERT_TRUE(e1.isEqualUndirected(e2));
    }
}

TEST_F(EdgeGTest, testWeightedEdgeBasics) {
    // print
    {
        std::stringstream ss;
        WeightedEdge e(4, 6, 7.5);
        ss << e;

        std::string s = ss.str();

        ASSERT_NE(s.find('4'), std::string::npos);
        ASSERT_NE(s.find('6'), std::string::npos);
        ASSERT_NE(s.find('7'), std::string::npos);
        ASSERT_NE(s.find('5'), std::string::npos);

        ASSERT_LE(s.find('7'), s.find('5'));
    }

    // undirected equality
    {
        WeightedEdge e1(3, 6, 12.0), e2(6, 3, 13.);

        ASSERT_FALSE(e1.isEqualDirected(e2));
        ASSERT_TRUE(e1.isEqualUndirected(e2));
        ASSERT_NE(e1, e2.flipped());
    }
}

} // namespace NetworKit
