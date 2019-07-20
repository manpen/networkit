/*
 * CurveballGlobalTradeGenerator.cpp
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include <cassert>
#include <algorithm>
#include <numeric>
#include <vector>

#include <networkit/randomization/CurveballGlobalTradeGenerator.hpp>

#include <networkit/auxiliary/Random.hpp>

namespace NetworKit {

CurveballGlobalTradeGenerator::value_type CurveballGlobalTradeGenerator::generate() const {
    value_type trades_out;
    trades_out.reserve(numGlobalTrades * numNodes / 2);

    std::vector<node> node_permutation(numNodes);
    std::iota(node_permutation.begin(), node_permutation.end(), 0);

    for (count run = 0; run < numGlobalTrades; run++) {
        // shuffling a shuffled node_permutation is okay, no need to reinitialize it
        std::shuffle(node_permutation.begin(), node_permutation.end(), Aux::Random::getURNG());

        // we do count the number of steps rather than checking iterator for end
        // since we do 2 increments per loop and the permutation vector can have odd length.
        // in this case we the last element (random node) is dropped.
        auto it = node_permutation.cbegin();
        for (count t_id = 0; t_id < numNodes / 2; t_id++) {
            assert(it != node_permutation.cend());

            const auto fst = *(it++);
            const auto snd = *(it++);
            trades_out.emplace_back(fst, snd);
        }
    }

    return trades_out;
}

} // namespace NetworKit
