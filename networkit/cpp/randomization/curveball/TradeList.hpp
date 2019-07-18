/*
 * TradeList.hpp
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */
#ifndef RANDOMIZATION_CURVEBALL_TRADE_LIST_H
#define RANDOMIZATION_CURVEBALL_TRADE_LIST_H

#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace CurveballDetails {

// Global Definitions
using trade_descriptor = std::pair<node, node>;
using tradeid = node;
using trade_vector = std::vector<trade_descriptor>;

static constexpr tradeid TRADELIST_END = std::numeric_limits<tradeid>::max();

class TradeList {
public:
    using edge_vector = std::vector<std::pair<node, node>>;
    using offset_vector = std::vector<tradeid>;
    using tradeid_vector = std::vector<tradeid>;
    using trade = trade_descriptor;
    using tradeid_it = tradeid_vector::const_iterator;

    TradeList(const node num_nodes) :
        numNodes(num_nodes)
    {}

    // Initialize method
    void initialize(const trade_vector &trades) {
        tradeList.clear();
        tradeList.resize(2 * trades.size() + numNodes);
        offsets.clear();
        offsets.resize(numNodes + 1);

        assert(numNodes > 0);
        assert(trades.size() > 0);

        std::vector<tradeid> trade_count(numNodes);

        // Push occurrences
        for (const auto trade : trades) {
            assert(trade.first < numNodes);
            assert(trade.second < numNodes);

            trade_count[trade.first]++;
            trade_count[trade.second]++;
        }

        // add missing +1 for sentinel
        trade_count[0]++;
        std::partial_sum(trade_count.cbegin(), trade_count.cend(),
                         offsets.begin() + 1,
                         [&](const tradeid a, const tradeid b) { return a + b + 1; });
        // add dummy
        offsets[numNodes] = 2 * trades.size() + numNodes - 1;

        // set sentinels
        for (node node = 1; node < numNodes; node++) {
            tradeList[offsets[node] - 1] = TRADELIST_END;
        }
        // set last entry as sentinel
        tradeList.back() = TRADELIST_END;

        std::fill(trade_count.begin(), trade_count.end(), 0);
        {
            tradeid trade_id = 0;
            for (const auto trade : trades) {
                auto updateNode = [&](const node u) {
                    const node pos = offsets[u] + trade_count[u];
                    tradeList[pos] = trade_id;
                    trade_count[u]++;
                };

                updateNode(trade.first);
                updateNode(trade.second);
                trade_id++;
            }
        }
    }

    // No Copy Constructor
    TradeList(const TradeList &) = delete;

    tradeid_it getTrades(const node nodeid) const {
        assert(nodeid < numNodes);

        return tradeList.begin() + offsets[nodeid];
    }

    void incrementOffset(const node nodeid) {
        assert(nodeid < numNodes);
        assert(1 <= offsets[nodeid + 1] - offsets[nodeid]);

        offsets[nodeid]++;
    }

    node numberOfNodes() const { return numNodes; }

protected:
    tradeid_vector tradeList;
    offset_vector offsets;
    const node numNodes;

};

} // namespace CurveballDetails
} // namespace NetworKit

#endif // RANDOMIZATION_CURVEBALL_TRADE_LIST_H
