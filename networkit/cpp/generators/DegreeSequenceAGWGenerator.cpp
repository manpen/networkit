#include <networkit/generators/DegreeSequenceAGWGenerator.hpp>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>

#include <algorithm>
#include <random>
#include <numeric>

namespace NetworKit {

#define ACCEPT return true
#define REJECT return false


DegreeSequenceAGWGenerator::DegreeSequenceAGWGenerator(const std::vector<count> &sequence) :
    StaticDegreeSequenceGenerator(sequence)
{}

Graph DegreeSequenceAGWGenerator::generate() {
    if (!num_nodes)
        return Graph(0);

    assert(isRealizable());
    initializeStats();

    while(true) {
        if (!sampleInitialMultiGraph()) continue;
        assert(edge_set.size() == num_edges - double_edges.size());

        if (!removeLoops()) continue;
        assert(nodes_with_loop.empty());

        if (!removeDoubleEdges()) continue;
        assert(double_edges.empty());
        assert(edge_set.size() == num_edges);

        break;
    }

    // free some memory early as it might be need to build the result graph
    node_stubs      = std::vector<node>();
    simple_degrees  = std::vector<count>();
    nodes_with_loop = std::vector<node>();
    double_edges    = std::vector<Edge>();

    return graphFromEdgeList();
}

void DegreeSequenceAGWGenerator::initializeStats() {
    // compute stats
    num_nodes     = seq.size();
    sum_degrees   = std::accumulate(begin(seq), end(seq), 0llu);
    num_edges     = sum_degrees/2;
    max_degree    = num_nodes ? *std::max_element(begin(seq), end(seq)) : 0;


    num_two_paths = std::accumulate(begin(seq), end(seq), 0llu, [] (count x) {two_path_from_degree(x);});

    num_max_loops = 22 * max_degree * max_degree * max_degree < num_two_paths ? num_two_paths / num_edges : 0;
    num_max_double_edges = num_max_loops * num_max_loops;

    // some plausibility checks on the distribution
    if (!num_nodes)
        return;

    if (max_degree >= num_nodes)
        throw std::runtime_error("Max degree matches or exceeds number of nodes");

    if (sum_degrees % 2)
        throw std::runtime_error("Sum of degrees is odd");

    // prepare sampling of configuration model: place d_i many stubs of node i into vector
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
}

bool DegreeSequenceAGWGenerator::sampleInitialMultiGraph() {
    std::random_shuffle(begin(node_stubs), end(node_stubs), Aux::Random::getURNG());

    // copy edges into edge vector
    edges.clear();
    edges.reserve(num_edges);
    nodes_with_loop.clear();
    for(auto it = begin(node_stubs); it != end(node_stubs); it += 2) {
        edges.emplace_back(*it, *(it + 1), true); // emplace sorted edge (i.e. first entry is smaller)
        if (edges.back().isLoop()) {
            nodes_with_loop.push_back(edges.back().u);

            if (nodes_with_loop.size() > num_max_loops)
                // too many self-loops
                REJECT;
        }
    }

    // TODO: Here we can use linear time integer sort
    Aux::Parallel::sort(begin(edges), end(edges));

    // scan through (sorted!) edge list to count double edges. if we find
    // a triple edge or double loop we reject and sample a new multigraph
    double_edges.clear();
    for (auto prev = begin(edges), cur = begin(edges) + 1; cur != end(edges); ++prev, ++cur) {
        if (*cur != *prev) {
            continue;
        }

        if (cur->isLoop()) {
            // we do not allow double self-loops
            REJECT;


        } else if (prev != begin(edges)) {
            double_edges.push_back(*cur);

            if (*std::prev(prev) == *cur || double_edges.size() > num_max_double_edges)
                // too many double edges
                REJECT;
        }
    }

    if (double_edges.empty() && nodes_with_loop.empty())
        // graph already valid
        ACCEPT;

    // build set
    edge_set.clear();
    edge_set.bulk_load(begin(edges), end(edges));

    assert(std::is_sorted(begin(nodes_with_loop), end(nodes_with_loop)));

    // compute number of two path
    {
        auto simple_degrees = seq;
        for(auto e : double_edges) {
            assert(simple_degrees[e.u] >= 2 & simple_degrees[e.v] >= 2);

            simple_degrees[e.u] -= 2;
            simple_degrees[e.v] -= 2;
        }

        // don't care about underflows, since if (x-1) underflows x will be zero and so will be x(x-1) ...
        // cast only to avoid warnings
        num_simple_twopaths = 0;
        auto loop_it = begin(nodes_with_loop);
        auto next_loop = (loop_it == end(nodes_with_loop) ? num_nodes : *loop_it);

        for(node u = 0; u < num_nodes; ++u) {
            if (u == next_loop) {
                // skip loop
                ++loop_it;
                next_loop = (loop_it == end(nodes_with_loop) ? num_nodes : *loop_it);
                continue;
            }

            num_simple_twopaths += two_path_from_degree(simple_degrees[u]);
        }

        assert(num_simple_twopaths <= num_two_paths);
    }

    ACCEPT;
}

bool DegreeSequenceAGWGenerator::removeLoops() {
    auto& gen = Aux::Random::getURNG();

    std::uniform_int_distribution<size_t> edge_distr{0, num_edges-1};

    while(!nodes_with_loop.empty()) {
        // draw random positions
        const auto loop_idx = std::uniform_int_distribution<size_t>{0, nodes_with_loop.size()-1}(gen);
        const auto edge1_idx = edge_distr(gen);
        const auto edge2_idx = edge_distr(gen);

        // load data
        const auto loop = nodes_with_loop[loop_idx];
        const auto edge1 = edges[edge1_idx];
        const auto edge2 = edges[edge2_idx];

        // reject if we did not obtain 5 different nodes
        if (edge1.isIncident(loop))  REJECT;
        if (edge1.isIncident(loop))  REJECT;
        if (edge1.isIncident(edge2)) REJECT;

        std::array<Edge, 3> new_edges{
            Edge{loop,    edge1.u, true},
            Edge{loop,    edge2.u, true},
            Edge{edge1.v, edge2.v, true}
        };

        if (std::any_of(begin(new_edges), end(new_edges), [this] (const Edge e) {return edge_set.count(e);})) {
            // if any new edge already exists: skip
            REJECT;
        }

        // we reached this point with probability f_l(G) / f_l(V_2). We completed the f-rejection.
        num_simple_twopaths += two_path_from_degree(simple_degrees[loop]);
        assert(num_simple_twopaths <= num_two_paths);

    }

    ACCEPT;
}

bool DegreeSequenceAGWGenerator::removeDoubleEdges() {
    auto& gen = Aux::Random::getURNG();

    std::uniform_int_distribution<size_t> edge_distr{0, num_edges-1};

    while(!double_edges.empty()) {

    }

    ACCEPT;
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