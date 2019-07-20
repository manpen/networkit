/*
 * CurveballUniformTradeGenerator.cpp
 *
 *  Created on: Jul 11, 2017
 *  	Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <cassert>
#include <algorithm>
#include <vector>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/randomization/CurveballUniformTradeGenerator.hpp>


namespace NetworKit {

CurveballUniformTradeGenerator::value_type CurveballUniformTradeGenerator::generate() const {
	value_type trades_out(numTrades);

	#pragma omp parallel if (numTrades > 10000)
	{
		auto& gen = Aux::Random::getURNG();

		#pragma omp for
		for (omp_index t_id = 0; t_id < static_cast<omp_index>(numTrades); ++t_id) {
			trades_out[t_id] = randomTrade(gen);
		}
	}

	return trades_out;
}

} // namespace NetworKit
