/*
 * CurveballUniformTradeGenerator.cpp
 *
 *  Created on: Jul 11, 2017
 *  	Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */

#include <cassert>
#include <algorithm>
#include <vector>

#include <tlx/define/likely.hpp>

#include <networkit/auxiliary/Random.hpp>

#include <networkit/randomization/CurveballUniformTradeGenerator.hpp>


namespace NetworKit {

CurveballUniformTradeGenerator::value_type CurveballUniformTradeGenerator::generate() const {
	value_type trades_out(numTrades);

	#pragma omp parallel if (numTrades > 10000)
	{
		const bool sampleFromVector = !allowedIds.empty();
		auto& gen = Aux::Random::getURNG();
		auto distr = std::uniform_int_distribution<node>(idLowerBound, idUpperBound);

		auto randomNode = [&] {
			return sampleFromVector ? allowedIds[distr(gen)] : distr(gen);
		};

		#pragma omp for
		for (omp_index t_id = 0; t_id < static_cast<omp_index>(numTrades); ++t_id) {
			const node fst = randomNode();
			node snd;
			do {
				snd = randomNode();
			} while (TLX_UNLIKELY(fst == snd));

			trades_out[t_id] = {fst, snd};
		}
	}

	return trades_out;
}

} // namespace NetworKit
