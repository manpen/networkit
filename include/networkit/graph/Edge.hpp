/*
 * Edge.hpp
 *
 *  Created on: 13.09.2019
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef NETWORKIT_GRAPH_EDGE_HPP
#define NETWORKIT_GRAPH_EDGE_HPP

#include <algorithm>
#include <utility>
#include <tuple>
#include <string>
#include <ostream>

#include <networkit/Globals.hpp>

namespace NetworKit {

template <typename T>
class EdgeLogic {
public:
// manipulate
    ///! Ensures that first <= second.
    T& sort() noexcept {
        if (dref().first > dref().second) flip();
        return dref();
    }

    ///! Returns a copy of this edge but with first <= second.
    T sorted() const noexcept {
        return T{dref(), true};
    }

    ///! Flips head and tail node in place
    T& flip() noexcept {
        std::swap(dref().first, dref().second);
        return dref();
    }

    ///! Returns a new edge with head and tail flipped
    T flipped() const noexcept {
        return T{dref()}.flip();
    }

// checks
    ///! Test whether this edge is a loop, i.e. connected to only one nde
    bool isLoop() const noexcept {
        return dref().first == dref().second;
    }

    ///! Test whether this edge is incident to node x.
    bool isIncident(const node x) const noexcept {
        return dref().first == x || dref().second == x;
    }

    ///! Test whether this edge is incident to edge e, i.e. share at least one node
    bool isIncident(const T e) const noexcept {
        return isIncident(e.first) || isIncident(e.second);
    }

    /// Test whether both endpoints match respecting their order.
    /// This methods only considers topology and ignores additional attributes such as weights.
    bool isEqualDirected(const T& o) const noexcept {
        return std::tie(dref().first, dref().second) == std::tie(o.first, o.second);
    }

    /// Test whether both endpoints match ignoring their order
    /// This methods only considers topology and ignores additional attributes such as weights.
    bool isEqualUndirected(const T& o) const noexcept {
        return isEqualDirected(o) || isEqualDirected(o.flipped());
    }

private:
    T& dref() {
        return *static_cast<T*>(this);
    }

    const T& dref() const {
        return *static_cast<const T*>(this);
}
};

class Edge : public EdgeLogic<Edge> {
public:
    // we keep the names first and second for compatibility with std::pair<node, node>
    node first;  ///< Tail node
    node second; ///< Head node

    Edge() {}

    Edge(node first, node second, bool sorted = false) noexcept
        : first(first), second(second)
    {
        if (sorted) sort();
    }

    Edge(const Edge& e, bool sorted = false) noexcept :
        first(e.first), second(e.second)
    {
        if (sorted) sort();
    }

    Edge(Edge&& e, bool sorted = false) noexcept :
        first(e.first), second(e.second)
    {
        if (sorted) sort();
    }

    Edge &operator=(const Edge &) noexcept = default;
    Edge &operator=(Edge &&) noexcept = default;

// compatibility to std::pair and std::tuple
    explicit Edge(std::tuple<node, node> pair, bool sorted = false) noexcept :
        first(std::get<0>(pair)), second(std::get<1>(pair))
    {
        if (sorted) sort();
    }

    explicit Edge(std::pair<node, node> pair, bool sorted = false) noexcept :
        first(pair.first), second(pair.second)
    {
        if (sorted) sort();
    }

// compare
    /// Lexicographical compare of (first, second) == (o.first, o.second)
    bool operator==(const Edge &o) const noexcept { return asPair() == o.asPair(); }
    
    /// Lexicographical compare of (first, second) != (o.first, o.second)
    bool operator!=(const Edge &o) const noexcept { return asPair() != o.asPair(); }
    
    /// Lexicographical compare of (first, second) <= (o.first, o.second)
    bool operator<=(const Edge &o) const noexcept { return asPair() <= o.asPair(); }
    
    /// Lexicographical compare of (first, second) >= (o.first, o.second)
    bool operator>=(const Edge &o) const noexcept { return asPair() >= o.asPair(); }
    
    /// Lexicographical compare of (first, second) <  (o.first, o.second)
    bool operator< (const Edge &o) const noexcept { return asPair() <  o.asPair(); }
    
    /// Lexicographical compare of (first, second) >  (o.first, o.second)
    bool operator> (const Edge &o) const noexcept { return asPair() >  o.asPair(); }

// cast to tuple
    /**
     * Implicit cast to a tuple of references. Useful to tie to an Edge:
     *
     * Edge e(1,2);
     * node u, v;
     * std::tie(u, v) = e;
     */
    operator std::tuple<const node&, const node&>() const noexcept {return std::tuple<const node&, const node&>{first, second};}
    operator std::tuple<      node&,       node&>()       noexcept {return std::tuple<      node&,       node&>{first, second};}
    operator std::tuple<      node ,       node >() const noexcept {return std::tuple<      node ,       node >{first, second};}

    operator std::pair<const node&, const node&>() const noexcept {return std::pair<const node&, const node&>{first, second};}
    operator std::pair<      node&,       node&>()       noexcept {return std::pair<      node&,       node&>{first, second};}
    operator std::pair<      node ,       node >() const noexcept {return std::pair<      node ,       node >{first, second};}


// print
    friend std::ostream& operator<<(std::ostream& os, const Edge& e) {
        const std::string s {"(" + std::to_string(e.first) + ", " + std::to_string(e.second) + ")"};
        return os << s;
    }

private:
    std::pair<node, node> asPair() const noexcept  {
        return {first, second};
    }

};

class WeightedEdge : public EdgeLogic<WeightedEdge> {
public:
    // we keep the names first and second for compatibility with std::pair<node, node>
    node first;  ///< Tail node
    node second; ///< Head node
    edgeweight weight;

    WeightedEdge() {} 

    WeightedEdge(node first, node second, edgeweight w, bool sorted = false) noexcept :
        first(first), second(second), weight(w)
    {
        if (sorted) sort();
    }

    WeightedEdge(const WeightedEdge& edge, bool sorted = false) noexcept :
        first(edge.first),
        second(edge.second),
        weight(edge.weight)
    {
        if (sorted) sort();
    }

    WeightedEdge(WeightedEdge&& edge, bool sorted = false) noexcept :
        first(edge.first),
        second(edge.second),
        weight(edge.weight)
    {
        if (sorted) sort();
    }

    WeightedEdge &operator=(const WeightedEdge &) noexcept = default;
    WeightedEdge &operator=(WeightedEdge &&) noexcept = default;

// compatibilty to std::pair and std::tuple
    explicit WeightedEdge(std::tuple<node, node, edgeweight> edge, bool sorted = false) noexcept :
        first(std::get<0>(edge)),
        second(std::get<1>(edge)),
        weight(std::get<2>(edge))
    {
        if (sorted) sort();
    }

// compare
    /// Lexicographical compare of (weight, first, second) == (o.weight, o.first, o.second)
    bool operator==(const WeightedEdge &o) const noexcept { return asTuple() == o.asTuple(); }

    /// Lexicographical compare of (weight, first, second) != (o.weight, o.first, o.second)
    bool operator!=(const WeightedEdge &o) const noexcept { return asTuple() != o.asTuple(); }

    /// Lexicographical compare of (weight, first, second) <= (o.weight, o.first, o.second)
    bool operator<=(const WeightedEdge &o) const noexcept { return asTuple() <= o.asTuple(); }

    /// Lexicographical compare of (weight, first, second) >= (o.weight, o.first, o.second)
    bool operator>=(const WeightedEdge &o) const noexcept { return asTuple() >= o.asTuple(); }

    /// Lexicographical compare of (weight, first, second) <  (o.weight, o.first, o.second)
    bool operator< (const WeightedEdge &o) const noexcept { return asTuple() <  o.asTuple(); }

    /// Lexicographical compare of (weight, first, second) >  (o.weight, o.first, o.second)
    bool operator> (const WeightedEdge &o) const noexcept { return asTuple() >  o.asTuple(); }

// casts
    operator   std::tuple<const node&, const node&, const edgeweight&>() const noexcept {
        return std::tuple<const node&, const node&, const edgeweight&>{first, second, weight};}
    operator   std::tuple<      node&,       node&,       edgeweight&>()       noexcept {
        return std::tuple<      node&,       node&,       edgeweight&>{first, second, weight};}
    operator   std::tuple<      node,        node,        edgeweight >() const noexcept {
        return std::tuple<      node,        node,        edgeweight >{first, second, weight};}

// print
    friend std::ostream& operator<<(std::ostream& os, const WeightedEdge& e) {
        const std::string s {"(" + std::to_string(e.first) + ", " + std::to_string(e.second) + ", " + std::to_string(e.weight) + ")"};
        return os << s;
    }

private:
    std::tuple<edgeweight, node, node> asTuple() const noexcept {
        return {weight, first, second};
    }
};

} // ! namespace NetworKit

#endif // NETWORKIT_GRAPH_EDGE_HPP
