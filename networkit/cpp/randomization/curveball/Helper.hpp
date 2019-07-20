/*
 * CurveballImpl.h
 *
 * Author: Manuel Penschuck <networkit@manuel.jetzt>
 */
#ifndef RANDOMIZATION_CURVEBALL_HELPER_H
#define RANDOMIZATION_CURVEBALL_HELPER_H

#include <algorithm>
#include <cassert>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {
namespace CurveballDetails {

template<typename UIterator, typename VIterator>
void computeCommonDisjointNeighbour(
    UIterator ubegin, UIterator uend,
    VIterator vbegin, VIterator vend,
    std::vector<node> &common_neighbours,
    std::vector<node> &disjoint_neighbours)
{
    if (std::distance(ubegin, uend) > std::distance(vbegin, vend))
        return computeCommonDisjointNeighbour(vbegin, vend, ubegin, uend, common_neighbours, disjoint_neighbours);

    constexpr node BIT = node(1) << (sizeof(node) * 8 - 1);
    constexpr node MASK = ~BIT;

    common_neighbours.clear();
    disjoint_neighbours.clear();

    std::sort(ubegin, uend);

    size_t remaining_hits = std::distance(ubegin, uend) / 2;

    #ifndef NDEBUG
    const size_t initial_size = std::distance(ubegin, uend) + std::distance(vbegin, vend);
    #endif

    for(; vbegin != vend; ++vbegin) {
        const auto nv = *vbegin;
        const auto u_it = std::lower_bound(ubegin, uend, nv, [MASK] (const node u, const node v) {return (u&MASK) < v;});

        if (u_it != uend && *u_it == nv) {
            common_neighbours.push_back(nv);
            *u_it |= BIT;

            if (!--remaining_hits)
            {
                uend = std::remove_if(ubegin, uend,
                                     [BIT] (const node u) {return u & BIT;});

                remaining_hits = std::distance(ubegin, uend) / 2;
                if (remaining_hits < 8)
                    remaining_hits = std::distance(ubegin, uend);
            }

        } else {
            disjoint_neighbours.push_back(nv);
        }
    }

    uend = std::remove_if(ubegin, uend, [BIT] (const node u) {return u & BIT;});
    disjoint_neighbours.insert(disjoint_neighbours.end(), ubegin, uend);

    assert(2*common_neighbours.size() + disjoint_neighbours.size() == initial_size);
    #ifndef NDEBUG
    for(auto x : common_neighbours)
        assert(std::find(disjoint_neighbours.cbegin(), disjoint_neighbours.cend(), x) == disjoint_neighbours.cend());
    #endif
}

} // namespace CurveballDetails
} // namespace NetworKit

#endif // RANDOMIZATION_CURVEBALL_HELPER_H
