/*
 * AdjacencyList.hpp
 *
 * Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */
#ifndef RANDOMIZATION_CURVEBALL_ADJCENCY_LIST_H
#define RANDOMIZATION_CURVEBALL_ADJCENCY_LIST_H

#include <algorithm>
#include <cassert>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace CurveballDetails {

class AdjacencyList {
public:
    using degree_vector = std::vector<count>;
    using neighbour_vector = std::vector<node>;
    using pos_vector = std::vector<edgeid>;
    using pos_it = pos_vector::iterator;
    using neighbour_it = neighbour_vector::iterator;
    using cneighbour_it = neighbour_vector::const_iterator;
    using nodepair_vector = std::vector <std::pair<node, node>>;

    AdjacencyList() = delete;

    AdjacencyList(const Graph &G) :
        neighbours(2 * G.numberOfEdges() + G.numberOfNodes() + 1),
        offsets(G.numberOfNodes()),
        begins(G.numberOfNodes() + 1)
    {
        degreeCount = 2 * G.numberOfEdges();

        count sum = 0;

        G.forNodes([&](node u) {
            begins[u] = sum;

            sum += G.degree(u);
            neighbours[sum] = LISTROW_END;

            // shift after Sentinel
            sum += 1;
        });

        neighbours[sum] = LISTROW_END;
        begins[G.numberOfNodes()] = sum;

        assert(sum == neighbours.size() - 1);
    }

    // no copy
    AdjacencyList(const AdjacencyList &) = delete;
    AdjacencyList &operator=(const AdjacencyList &) = delete;

    // movable
    AdjacencyList(AdjacencyList &&) = default;
    AdjacencyList &operator=(AdjacencyList &&) = default;

    void restructure() {
        std::fill(offsets.begin(), offsets.end(), 0);
    }

    template<typename CB>
    void forEdges(CB &&callback) const {
        for (node nodeid = 0; nodeid < static_cast<node>(offsets.size()); nodeid++) {
            for (auto it = cbegin(nodeid); it != cend(nodeid); it++) {
                callback(nodeid, *it);
            }
        }
    }

    nodepair_vector getEdges() const {
        nodepair_vector edges;
        edges.reserve(degreeCount);
        forEdges([&edges](node u, node v) {
            edges.emplace_back(u, v);
        });

        return edges;
    }

    void insertNeighbour(const node node_id, const node neighbour) {
        auto pos = begin(node_id) + offsets[node_id];

        assert(*pos != LISTROW_END);

        *pos = neighbour;

        offsets[node_id]++;
    }

    node numberOfNodes() const {
        return static_cast<node>(offsets.size());
    }

    node numberOfEdges() const {
        return static_cast<edgeid>(degreeCount);
    }

    void resetRow(const node node_id) {
        assert(node_id < static_cast<node>(offsets.size()));
        offsets[node_id] = 0;
    }

    count degreeAt(node node_id) const {
        assert(node_id < static_cast<node>(offsets.size()));
        return begins[node_id + 1] - begins[node_id] - 1;
    }

    neighbour_it begin(const node node_id) {
        return neighbours.begin() + begins[node_id];
    }

    neighbour_it end(const node node_id) {
        return begin(node_id) + offsets[node_id];
    }

    cneighbour_it cbegin(const node node_id) const {
        return neighbours.cbegin() + begins[node_id];
    }

    cneighbour_it cend(const node node_id) const {
        return cbegin(node_id) + offsets[node_id];
    }

protected:
    static constexpr count LISTROW_END = std::numeric_limits<count>::max();

    neighbour_vector neighbours;
    degree_vector offsets;
    pos_vector begins;
    edgeid degreeCount;

};

} // namespace CurveballDetails
} // namespace NetworKit


#endif // RANDOMIZATION_CURVEBALL_ADJCENCY_LIST_H
