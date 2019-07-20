/*
 * CurveballUniformTradeGenerator.h
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */
#ifndef RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_H
#define RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_H

#include <algorithm>
#include <cassert>
#include <vector>
#include <random>
#include <utility>

#include <tlx/define/likely.hpp>

#include <networkit/Globals.hpp>

namespace NetworKit {

/**
 * Generate a sequence of uniform Curveball trades
 *
 * A Curveball trade is expressed as a pair of nodes. Each node is drawn
 * independently and uniformly at random. Self-trades of form (x, x) are
 * forbidden.
 */
class CurveballUniformTradeGenerator {
public:
    using trade = std::pair<node, node>;
    using value_type = std::vector<trade>;

protected:
	count numTrades;

	node idLowerBound;
	node idUpperBound;

	std::vector<node> allowedIds;

public:
	///! Trade @a trade_num many trades of nodes from the id range of [@a idLowerBound, @a idUpperBound).
	CurveballUniformTradeGenerator(const count trade_num, const node idLowerBound, const node idUpperBound)
		: numTrades(trade_num), idLowerBound(idLowerBound), idUpperBound(idUpperBound - 1)
	{
		if (idLowerBound + 1 >= idUpperBound)
			throw std::runtime_error("idUpperBound has to strictly greater than idLowerBound with a distance of at least 2");
	}

	///! Trade @a trade_num many trades of nodes from the id range of [0, @a idUpperBound).
	CurveballUniformTradeGenerator(const count trade_num, const node idUpperBound)
		: CurveballUniformTradeGenerator(trade_num, 0, idUpperBound)
	{}

	/***
	 * Trade @a trade_num many trades of nodes the set @a allowedIds. Entries are sampled
	 * uniformly at random, so if a node appears multiple times in @a allowedIds,
	 * it has a higher chance to be sampled.
	 *
	 * @warning @a allowedIds has to contain at least two different ids
	 */
	CurveballUniformTradeGenerator(const count trade_num, std::vector<node> allowedNodes)
		: numTrades(trade_num), idLowerBound(0), idUpperBound(allowedNodes.size() - 1), allowedIds(std::move(allowedNodes))
	{
		if (allowedIds.size() < 2)
			throw std::runtime_error("Need at least two allowed ids");

		// make sure that there are at least two different entries in allowedIds
		assert(std::find_if(allowedIds.cbegin(), allowedIds.cend(), [&] (node x) {return x != allowedIds.front();}) != end(allowedIds));
	}

	///! Generate and return a trade sequence. May be called multiple times.
	value_type generate() const;

	///! Returns two different random nodes
	trade randomTrade(std::mt19937_64& gen) const {
	// Helper to sample random node
		auto distr = std::uniform_int_distribution<node>(idLowerBound, idUpperBound);
		auto randomNode = [&] {
			return allowedIds.empty() ? distr(gen) : allowedIds[distr(gen)];
		};

	 // Sample pair of different nodes
		const node fst = randomNode();
		node snd;
		do {
			snd = randomNode();
		} while (TLX_UNLIKELY(fst == snd));

		return {fst, snd};
	}

	count numberOfTrades() const {
		return numTrades;
	}
};

}

#endif // RANDOMIZATION_CURVEBALL_UNIFORM_TRADE_GENERATOR_H
