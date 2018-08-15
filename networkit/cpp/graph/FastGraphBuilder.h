#include <omp.h>
#include <atomic>
#include <iterator>
#include <memory>
#include <mutex>
#include <vector>

#include "Graph.h"
#include "../Globals.h"

namespace NetworKit {
namespace FastGraphBuilderDetails {
    template<bool WithWeights> class Message; // forward declaration (see at the end of this file)

    template <typename T, size_t BlockSize>
    class Block {
    public:
        static constexpr size_t block_size = BlockSize;

        Block() {bufferEnd = &buffer[0];}

        T* begin() {return &buffer[0];}
        T* end() {return bufferEnd;}
        const T* begin() const {return &buffer[0];}
        const T* end() const {return bufferEnd;}

        size_t size() const {return std::distance(&buffer[0], const_cast<const T*>(bufferEnd));}
        bool empty() const {return buffer == bufferEnd;}
        bool full() const {return size() == block_size;}

        template <typename ... Args>
        void emplace_back(Args... args) {
            new (bufferEnd) T(std::forward<Args>(args)...);
            bufferEnd++;
        }

    private:
        T buffer[block_size];
        T* bufferEnd;
    };
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
    {
        for(size_t tid=0; tid != insertionBuffers.size(); tid++)
            insertionBuffers[tid].reset(new block_type);
    }

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
        assert(insertionBuffers[tid]);
        assert(u < n);
        assert(v < n);

        auto& buffer = insertionBuffers[tid];
        buffer->emplace_back(u, v, false, weight);

        if (buffer->full()) {
            {
                std::unique_lock<std::mutex> lock(filledBlocksMutex);
                filledBlocks.emplace_back(std::move(buffer));
            }

            buffer.reset(new block_type);
        }
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
    using block_type = FastGraphBuilderDetails::Block<message_type, 1024*1024>;
    std::vector< std::unique_ptr<block_type> > insertionBuffers;

    std::mutex filledBlocksMutex;
    std::list< std::unique_ptr<block_type> > filledBlocks;

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