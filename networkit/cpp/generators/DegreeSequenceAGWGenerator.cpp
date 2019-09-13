#include <networkit/generators/DegreeSequenceAGWGenerator.hpp>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>

#include <algorithm>
#include <random>
#include <numeric>

namespace NetworKit {


DegreeSequenceAGWGenerator::DegreeSequenceAGWGenerator(const std::vector<count> &sequence) :
    StaticDegreeSequenceGenerator(sequence),
    num_nodes(sequence.size()),
    sum_degrees(std::accumulate(begin(sequence), end(sequence), 0llu)),
    num_edges(sum_degrees/2),
    max_degree(num_nodes ? *std::max_element(begin(sequence), end(sequence)) : 0),
    num_two_paths(std::accumulate(begin(sequence), end(sequence), 0llu, [] (count x) {return x * (x-1);})),
    num_max_loops(22 * max_degree * max_degree * max_degree < num_two_paths ? num_two_paths / num_edges : 0),
    num_max_double_edges(num_max_loops * num_max_loops)
{
    if (!num_nodes)
        return;

    if (max_degree >= num_nodes)
        throw std::runtime_error("Max degree matches or exceeds number of nodes");

    if (sum_degrees % 2)
        throw std::runtime_error("Sum of degrees is odd");
}


Graph DegreeSequenceAGWGenerator::generate() {
    if (!num_nodes)
        return Graph(0);

    // place d_i many stubs of node i into vector
    node_stubs.resize(sum_degrees);
    {
        auto it = begin(node_stubs);
        node u = 0;
        for(const auto d : seq) {
            std::fill_n(it, d, u);
            it += d;
            ++u;
        }
        assert(it == end(node_stubs));
    }

    while(true) {
        sampleInitialMultiGraph(num_max_loops, num_max_double_edges);

    }

    return graphFromEdgeList();
}

bool DegreeSequenceAGWGenerator::sampleInitialMultiGraph(NetworKit::count allowed_num_loops,
                                                              NetworKit::count allowed_num_double_edges) {


    while(true) {
        std::random_shuffle(begin(node_stubs), end(node_stubs), Aux::Random::getURNG());

        // copy edges into edge vector
        edges.clear();
        edges.reserve(num_edges);
        nodes_with_loop.clear();
        for(auto it = begin(node_stubs); it != end(node_stubs); it += 2) {
            edges.emplace_back(*it, *(it + 1), true); // emplace sorted edge (i.e. first entry is smaller)
            if (edges.back().isLoop()) {
                nodes_with_loop.push_back(edges.back().u);
                if (nodes_with_loop.size() > allowed_num_loops)
                    break;
            }
        }

        // break early if there are too many self loops
        if (nodes_with_loop.size() > allowed_num_loops)
            continue;

        // TODO: Here we can use linear time integer sort
        Aux::Parallel::sort(begin(edges), end(edges));

        // scan through (sorted!) edge list to count double edges. if we find
        // a triple edge or double loop we reject and sample a new multigraph
        double_edges.clear();
        {
            bool reject = false;

            for (auto prev = begin(edges), cur = begin(edges) + 1; cur != end(edges); ++prev, ++cur) {
                if (*cur != *prev) {
                    continue;
                }

                if (cur->isLoop()) {
                    // we do not allow double self-loops
                    reject = true;
                    break;

                } else if (prev != begin(edges)) {
                    double_edges.push_back(*cur);

                    if (*std::prev(prev) == *cur || double_edges.size() > allowed_num_double_edges) {
                        reject = true;
                        break;
                    }
                }
            }

            if (reject)
                continue;
        }

        return;
    }
}

void DegreeSequenceAGWGenerator::removeLoops() {

}


Graph DegreeSequenceAGWGenerator::graphFromEdgeList() const {
    Graph G(num_nodes);

    for(node i = 0; i < num_nodes; ++i) {
        G.preallocateUndirected(i, seq[i]);
    }

    for(const Edge e : edges) {
        G.addPartialEdge(unsafe, e.u, e.v);
        G.addPartialEdge(unsafe, e.v, e.u);
    }

    G.setEdgeCount(unsafe, num_edges);

    return G;
}

}