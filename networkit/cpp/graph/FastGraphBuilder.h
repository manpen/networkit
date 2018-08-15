#include <omp.h>
#include <vector>
#include <deque>
#include <atomic>

#include "Graph.h"
#include "../Globals.h"

namespace NetworKit {
namespace FastGraphBuilderDetails {
    template<bool WithWeights> class Message; // forward declaration (see at the end of this file)
}


template <bool SupportWeight>
class FastGraphBuilder {
public:
    static constexpr bool supportWeights = SupportWeight;

    FastGraphBuilder() = delete;

    FastGraphBuilder(node n, bool weighted, bool directed, int numThreads = omp_get_max_threads())
        : insertionBuffers(numThreads),
          n(n),
          weighted(weighted),
          directed(directed)
    {}

    // not copyable
    FastGraphBuilder(const FastGraphBuilder&) = delete;
    FastGraphBuilder& operator=(const FastGraphBuilder&) = delete;

    // movable
    FastGraphBuilder(FastGraphBuilder&&) = default;
    FastGraphBuilder& operator=(FastGraphBuilder&&) = default;

    /**
     * Add an edge (u, v) -- and in case of undirected graphs also (v, u) --
     * to the graph. If the template parameter @a SupportWeight == false,
     * the weight parameter is silently ignored.
     *
     * @warning Do NOT call this function twice for the same edge; especially NOT
     * once per orientation of an undirected edge.
     *
     * @note While it's possible to call addEdge using omp_get_thread_num() as a
     * parameter, it's highly encouraged to cache its value for performance reasons.
     *
     * The edge is temporarily placed into the insertion buffer of thread tid.
     * @warning Ensure that no two threads add edges using the same tid to avoid
     * data races.
     */
    void addEdge(int tid, node u, node v, edgeweight weight = 1.0) {
        assert(tid < static_cast<int>(insertionBuffers.size()));
        assert(u < n);
        assert(v < n);

        insertionBuffers[tid].emplace_back(u, v, false, weight);
    }

    /**
     * Build and return the graph.
     *
     * @warning Calling this function invalidates the Graph Builders state.
     * Destroy it afterwards.
     */
    Graph toGraph();

private:
// Insertion Buffers
    using message_type = FastGraphBuilderDetails::Message<SupportWeight>;
    std::vector< std::deque<message_type> > insertionBuffers;
    std::vector< message_type > messages;

    void consolidateBuffers();

    node n;
    bool weighted;
    bool directed;
    std::atomic<count> selfLoops{0};
};

namespace FastGraphBuilderDetails {
class MessageBase {
public:
    MessageBase() {} // avoid initialisation when default constructing
    MessageBase(node tail, node head, bool flag) :
        tailNode(tail), headNodeFlag( (head << 1) | flag )
    {
        assert(tail < (1llu << (8*sizeof(node)-1)));
    }

    MessageBase(const MessageBase&) = default;

    node tail() const { return tailNode; }
    node head() const { return headNodeFlag >> 1; }
    bool flag() const { return headNodeFlag  & 1; }

private:
    // using a bit field causes errors in gcc (why?). Hence we use
    // explicit masking and shifting.
    node tailNode;
    node headNodeFlag;
};



template<>
class Message<false> : public MessageBase {
public:
    Message() {}
    Message(node tail, node head, bool flag, edgeweight = 1.0) : MessageBase(tail, head, flag) {}
    edgeweight weight() const { return 1.0; };
};

template<>
class Message<true> : public MessageBase {
public:
    Message() {}
    Message(node tail, node head, bool flag, edgeweight weight = 1.0) : MessageBase(tail, head, flag), eweight(weight) {}
    edgeweight weight() const { return eweight; };

private:
    edgeweight eweight;
};
}

}