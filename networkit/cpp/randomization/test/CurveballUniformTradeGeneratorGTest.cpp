/*
 * CurveballUniformTradeGeneratorGTest.h
 *
 *  Created on: 26.05.2018
 *      Author:  Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include <gtest/gtest.h>

#include <networkit/randomization/CurveballUniformTradeGenerator.hpp>

namespace NetworKit {

class CurveballUniformTradeGeneratorGTest : public ::testing::Test {};

TEST_F(CurveballUniformTradeGeneratorGTest, testGeneration) {
	CurveballUniformTradeGenerator gen(100, 10);
	auto trades = gen.generate();
	ASSERT_EQ(trades.size(), 100u);

	for (const auto trade : trades) {
		ASSERT_LE(trade.first, 9u);
		ASSERT_GE(trade.first, 0u);
		ASSERT_LE(trade.second, 9u);
		ASSERT_GE(trade.second, 0u);
	}

	ASSERT_TRUE(std::find_if(trades.cbegin(), trades.cend(), [] (std::pair<node, node> p) {return p.first == 0 || p.second == 0;}) != trades.cend());
	ASSERT_TRUE(std::find_if(trades.cbegin(), trades.cend(), [] (std::pair<node, node> p) {return p.first == 9 || p.second == 9;}) != trades.cend());
}

TEST_F(CurveballUniformTradeGeneratorGTest, testGenerationLB) {
	CurveballUniformTradeGenerator gen(105, 5, 10);
	auto trades = gen.generate();
	ASSERT_EQ(trades.size(), 105u);

	for (const auto trade : trades) {
		ASSERT_LE(trade.first, 9u);
		ASSERT_GE(trade.first, 5u);
		ASSERT_LE(trade.second, 9u);
		ASSERT_GE(trade.second, 5u);
	}

	ASSERT_TRUE(std::find_if(trades.cbegin(), trades.cend(), [] (std::pair<node, node> p) {return p.first == 5 || p.second == 5;}) != trades.cend());
	ASSERT_TRUE(std::find_if(trades.cbegin(), trades.cend(), [] (std::pair<node, node> p) {return p.first == 9 || p.second == 9;}) != trades.cend());
}

TEST_F(CurveballUniformTradeGeneratorGTest, testGenerationRestricted) {
	std::vector<node> allowed{5,6, 8,9}; // 7 is missing !

	CurveballUniformTradeGenerator gen(99, std::move(allowed));
	auto trades = gen.generate();

	ASSERT_EQ(trades.size(), 99u);

	for (const auto trade : trades) {
		ASSERT_LE(trade.first, 9u);
		ASSERT_GE(trade.first, 5u);
		ASSERT_NE(trade.first, 7u);
		ASSERT_LE(trade.second, 9u);
		ASSERT_GE(trade.second, 5u);
		ASSERT_NE(trade.second, 7u);
	}

	ASSERT_TRUE(std::find_if(trades.cbegin(), trades.cend(), [] (std::pair<node, node> p) {return p.first == 5 || p.second == 5;}) != trades.cend());
	ASSERT_TRUE(std::find_if(trades.cbegin(), trades.cend(), [] (std::pair<node, node> p) {return p.first == 9 || p.second == 9;}) != trades.cend());
}

}
