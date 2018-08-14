#include "FastGraphBuilder.h"

#include <algorithm>
#include <atomic>
#include "../auxiliary/IntSort.h"
#include "../auxiliary/Timer.h"

#include <parallel/algorithm>

namespace NetworKit {

template<bool SupportWeights>
void FastGraphBuilder<SupportWeights>::consolidateBuffers() {
    if (!messages.empty()) {
        #ifndef NDEBUG
        for(const auto& buf : insertionBuffers)
            assert(buf.empty());
        #endif

        return;
    }

    std::vector<size_t> begins;
    begins.reserve(insertionBuffers.size());

    size_t size = 0;
    for(const auto& buf : insertionBuffers) {
        begins.push_back(size);
        size += buf.size();
    }
    assert(size > 0);

    // copy partial insertion buffers into single large vector
    messages.resize(2*size);

    count num_selfloops = 0;
    #pragma omp parallel for reduction(+:num_selfloops)
    for(omp_index i=0; i<insertionBuffers.size(); i++) {
        if (insertionBuffers[i].empty())
            continue;

        auto writer = &messages[2*begins[i]];
        for(const auto& msg : insertionBuffers[i]) {
            *(writer++) = msg;
            //add a mirrored copy
            *(writer++) = message_type(msg.head(), msg.tail(), true, msg.weight());
            num_selfloops += (msg.head() == msg.tail());
        }
    }

    insertionBuffers.clear();
}

namespace FastGraphBuilderDetails {
    template <typename Callback, bool SupportWeights>
    void foreachNodeParallel(const std::vector<Message<SupportWeights> >& messages, Callback cb) {
        using citer = typename std::vector<Message<SupportWeights> >::const_iterator;

        #pragma omp parallel if (messages.size() > 1000)
        {
            const int threads = omp_get_num_threads();
            const int tid = omp_get_thread_num();

            const size_t chunk_size = (messages.size() + threads - 1) / threads; // ceildiv
            citer it                = messages.cbegin() + std::min<size_t>(messages.size(), chunk_size *  tid);
            const citer chunk_end   = messages.cbegin() + std::min<size_t>(messages.size(), chunk_size * (tid + 1));

            assert(tid != threads-1 || chunk_end == messages.cend());

            if (it != chunk_end) {
                auto next_node = [&messages](citer it) {
                    node u = it->tail();
                    for (++it; it != messages.cend() && it->tail() == u; ++it) {}
                    return it;
                };

                // begin of first thread is fixed, all other advance to their first "own" node,
                // i.e., if the last node of the previous chunk overlaps the boundary, it will be
                // processed by the previous PE
                if (tid && it->tail() == (it - 1)->tail())
                    it = next_node(it);

                while (it < chunk_end) {
                    citer node_end = next_node(it);
                    cb(it, node_end);
                    it = node_end;
                }
            }
        }

    }
}

template<bool SupportWeights>
Graph FastGraphBuilder<SupportWeights>::toGraph() {
    // compact insertion buffers
	consolidateBuffers();

	Aux::intsort(messages, [] (const message_type& m) {return m.tail();}, n);

    Graph G(n, weighted, directed);
    G.m = messages.size() / 2;
    G.storedNumberOfSelfLoops = selfLoops;

    using citer = typename std::vector<message_type>::const_iterator;
    if (directed) {
        FastGraphBuilderDetails::foreachNodeParallel<>(messages, [&](citer begin, citer end) {
            const size_t ideg = std::count_if(begin, end, [] (const message_type& m) {return m.flag();} );
            const size_t odeg = std::distance(begin, end) - ideg;
            const node u = begin->head();

            G.inDeg[u] = ideg;
            auto &inEdges = G.inEdges[u];
            auto &inWeights = G.inEdgeWeights[u];

            G.outDeg[u] = odeg;
            auto& outEdges = G.outEdges[u];
            auto& outWeights = G.outEdgeWeights[u];

            inEdges.reserve(ideg);
            outEdges.reserve(odeg);
            if (SupportWeights && weighted) {
                inWeights.reserve(ideg);
                outWeights.reserve(odeg);
            }

            for (auto it = begin; it != end; ++it) {
                if (it->flag()) {
                    inEdges.push_back(it->tail());
                    if (SupportWeights && weighted)
                        inWeights.push_back(it->weight());
                } else {
                    outEdges.push_back(it->tail());
                    if (SupportWeights && weighted)
                        outWeights.push_back(it->weight());
                }
            }
        });
    } else {
        using citer = typename std::vector<message_type>::const_iterator;
        FastGraphBuilderDetails::foreachNodeParallel<>(messages, [&] (citer begin, citer end) {
            const size_t deg = std::distance(begin, end);
            const node u = begin->tail();

            G.outDeg[u] = deg;

            auto& outEdges = G.outEdges[u];
            auto& outWeights = G.outEdgeWeights[u];

            outEdges.reserve(outEdges.size() + deg);
            if (SupportWeights && weighted)
                outWeights.reserve(deg);

            for(auto it=begin; it != end; ++it) {
                assert(it->tail() == u);
                outEdges.push_back(it->head());
                if (SupportWeights && weighted)
                    outWeights.push_back(it->weight());
            }
        });
    }

    return G;
}

template class FastGraphBuilder<false>;
template class FastGraphBuilder<true>;

}