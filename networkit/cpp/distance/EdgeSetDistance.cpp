#include <atomic>

#include <networkit/distance/EdgeSetDistance.hpp>
#include <tlx/define.hpp>

namespace NetworKit {

EdgeSetDistance::EdgeSetDistance(const Graph& g1, const Graph& g2)
    : Algorithm(), g1(g1), g2(g2)
{
    if  (g1.isDirected() != g2.isDirected())
        throw std::runtime_error("EdgeSetDistance requires that either both graphs are directed or both are undirected");
}

void EdgeSetDistance::run() {
    std::atomic<bool> node_mismatch{false};

    count num_edge_mismatch = 0;

    #pragma omp parallel reduction(+:num_edge_mismatch)
    {
        std::vector<node> n1;
        std::vector<node> n2;

        const auto upper = std::max<omp_index>(g1.upperNodeIdBound(), g2.upperNodeIdBound());
        const auto directed = g1.isDirected();

        #pragma omp for schedule(guided)
        for(omp_index u = 0; u < upper; ++u) {
            const bool exists1 = g1.hasNode(u);
            const bool exists2 = g2.hasNode(u);

            const auto deg = (exists1 ? g1.degree(u) : 0) + (exists2 ? g2.degree(u) : 0);
            if (!deg)
                continue;

            if (directed) {
                if (exists1 != exists2) {
                    num_edge_mismatch += deg;
                    continue;
                }

                n1.reserve(g1.degree(u));
                n2.reserve(g2.degree(u));

                g1.forNeighborsOf(u, [&](const node v) { n1.push_back(v); });
                g2.forNeighborsOf(u, [&](const node v) { n2.push_back(v); });

            } else {
                if (exists1 != exists2) {
                    if (exists1)
                        g1.forNeighborsOf(u, [&](const node v) { num_edge_mismatch += (u >= v); });
                    else
                        g2.forNeighborsOf(u, [&](const node v) { num_edge_mismatch += (u >= v); });

                    continue;
                }

                n1.reserve(g1.degree(u));
                n2.reserve(g2.degree(u));

                // in undirected graphs have to count only half of the edges
                g1.forNeighborsOf(u, [&](const node v) { if (u >= v) n1.push_back(v); });
                g2.forNeighborsOf(u, [&](const node v) { if (u >= v) n2.push_back(v); });
            }

            if (n1.empty()) {
                num_edge_mismatch += n2.size();
            } else if (n2.empty()) {
                num_edge_mismatch += n1.size();
            } else {
                std::sort(n1.begin(), n1.end());
                std::sort(n2.begin(), n2.end());

                // count (we could use std::set_intersection with an counting output iterator, but that's probabily more code)
                {
                    auto it1 = n1.cbegin();
                    auto it2 = n2.cbegin();

                    while(true) {
                        if (*it1 == *it2) {
                            ++it1; it2++;
                            if (it1 == n1.cend()) break;
                            if (it2 == n2.cend()) break;
                        } else {
                            ++num_edge_mismatch;

                            if (*it1 < *it2) {
                                if (++it1 == n1.cend()) break;
                            } else {
                                if (++it2 == n2.cend()) break;
                            }
                        }
                    }

                    num_edge_mismatch += std::distance(it1, n1.cend());
                    num_edge_mismatch += std::distance(it2, n2.cend());
                }
            }

            n1.clear();
            n2.clear();
        }
    };

    num_common_edges = (g1.numberOfEdges() + g2.numberOfEdges() - num_edge_mismatch) / 2;

    hasRun = true;
}

} // namespace NetworKit
