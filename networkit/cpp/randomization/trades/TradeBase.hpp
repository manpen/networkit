/*
 * Helper.hpp
 *
 * Author: Manuel Penschuck <networkit@manuel.jetzt>
 */
// networkit-format
#ifndef RANDOMIZATION_TRADES_HELPER_H
#define RANDOMIZATION_TRADES_HELPER_H

#include <range/v3/action/drop.hpp>

namespace NetworKit {
namespace CurveballDetails {
namespace Trades {

class TradeBase {
public:
    bool hasSharedEdge() const { return sharedEdgeType != SharedEdge::None; }
    bool hasSharedEdgeAtU() const { return sharedEdgeType == SharedEdge::AtU; }
    bool hasSharedEdgeAtV() const { return sharedEdgeType == SharedEdge::AtV; }

protected:
    template <typename RangeU, typename RangeV, typename CallbackCommon, typename CallbackDisjointU,
              typename CallbackDisjointV>
    void computeCommonDisjoint(RangeU u, RangeV v, CallbackCommon cb_common, CallbackDisjointU cb_u,
                               CallbackDisjointV cb_v) {

        auto uit = std::begin(u);
        auto vit = std::begin(v);

        if (!u.empty() && !v.empty()) {
            while (true) {
                if (*uit < *vit) {
                    cb_u(*uit);
                    if (++uit == u.end())
                        break;

                } else if (*vit < *uit) {
                    cb_v(*vit);
                    if (++vit == v.end())
                        break;

                } else {
                    cb_common(*uit);
                    ++uit;
                    ++vit;

                    if (uit == u.end() || vit == v.end())
                        break;
                }
            }
        }

        for (; uit != u.end(); ++uit)
            cb_u(*uit);
        for (; vit != v.end(); ++vit)
            cb_v(*vit);
    }

    enum class SharedEdge { None = 0, AtU = 1, AtV = 2 };

    template <typename RangeU, typename RangeV, typename CallbackCommon, typename CallbackDisjointU,
              typename CallbackDisjointV>
    void computeCommonDisjoint(node u, RangeU nu, node v, RangeV nv, CallbackCommon cb_common,
                               CallbackDisjointU cb_u, CallbackDisjointV cb_v) {

        sharedEdgeType = SharedEdge::None;
        computeCommonDisjoint(
            std::move(nu), std::move(nv), std::move(cb_common),
            [&](node x) {
                if (x == v) {
                    assert(sharedEdgeType
                           == SharedEdge::None); // may only have a single shared edge
                    sharedEdgeType = SharedEdge::AtU;
                } else {
                    cb_u(x);
                }
            },
            [&](node x) {
                if (x == u) {
                    assert(sharedEdgeType
                           == SharedEdge::None); // may only have a single shared edge
                    sharedEdgeType = SharedEdge::AtV;
                } else {
                    cb_v(x);
                }
            });
    }

    SharedEdge sharedEdgeType;
};

} // namespace Trades
} // namespace CurveballDetails
} // namespace NetworKit

#endif