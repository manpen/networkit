/*
 * PermutationTrade.hpp
 *
 * Author: Manuel Penschuck <networkit@manuel.jetzt>
 */
// networkit-format
#ifndef RANDOMIZATION_TRADES_RANDOM_PERMUTATION_H
#define RANDOMIZATION_TRADES_RANDOM_PERMUTATION_H

#include <algorithm>
#include <cassert>
#include <vector>

#include <range/v3/view.hpp>

#include <tlx/algorithm/random_bipartition_shuffle.hpp>
#include <tlx/define.hpp>
#include <tlx/die.hpp>

#include <networkit/Globals.hpp>

#include "TradeBase.hpp"

namespace NetworKit {
namespace CurveballDetails {
namespace Trades {

class RandomPerumation : public TradeBase {
public:
    explicit RandomPerumation(count maxDegree = 0) {
        commonNeighbors.reserve(maxDegree);
        disjointNeighbors.reserve(maxDegree);
    }

    template <typename RangeU, typename RangeV, typename Prng>
    void run(node u, RangeU neigh_u, node v, RangeV neigh_v, Prng &prng) {
        commonNeighbors.clear();
        disjointNeighbors.clear();

        computeCommonDisjoint(
            u, neigh_u, v, neigh_v,                          //
            [&](node x) { commonNeighbors.push_back(x); },   //
            [&](node x) { disjointNeighbors.push_back(x); }, //
            [&](node x) { disjointNeighbors.push_back(x); });

        disjointsOfU = ranges::distance(neigh_u) - commonNeighbors.size() - hasSharedEdgeAtU();
        assert(commonNeighbors.size() + ranges::distance(neighborsOfU()) + hasSharedEdgeAtU()
               == ranges::distance(neigh_u));
        assert(commonNeighbors.size() + ranges::distance(neighborsOfV()) + hasSharedEdgeAtV()
               == ranges::distance(neigh_v));

        tlx::random_bipartition_shuffle(disjointNeighbors.begin(), disjointNeighbors.end(),
                                        disjointsOfU, prng);
    }

    const auto &neighborsInCommon() const { return commonNeighbors; }
    auto neighborsOfU() const { return disjointNeighbors | ranges::views::take(disjointsOfU); }
    auto neighborsOfV() const { return disjointNeighbors | ranges::views::drop(disjointsOfU); }

private:
    std::vector<node> commonNeighbors;
    std::vector<node> disjointNeighbors;

    size_t disjointsOfU;
};

} // namespace Trades
} // namespace CurveballDetails
} // namespace NetworKit

#endif // ! RANDOMIZATION_TRADES_RANDOM_PERMUTATION_H