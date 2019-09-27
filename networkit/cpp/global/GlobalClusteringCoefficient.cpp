/*
 * GlobalClusteringCoefficient.cpp
 *
 *  Created on: 12.11.2013
 */

#include <networkit/global/GlobalClusteringCoefficient.hpp>
#include <networkit/auxiliary/Random.hpp>

#include <algorithm>

namespace NetworKit {

double GlobalClusteringCoefficient::approximate(const Graph& G, int k) {
  const count n = G.numberOfNodes();
  
  std::vector<int> w(n + 1);
  int sum = 0;
  for(node i = 0; i < n; i++) {
    w[i] = sum;
    sum += (G.degree(i) * (G.degree(i) - 1)) / 2;
  }
  w[n] = sum;

  int l = 0;
  for(count i = 0; i < k; i++) {
    int r2 = Aux::Random::index(w[n]);
    node r = static_cast<node>(std::distance(w.cbegin(),
        std::lower_bound(w.cbegin(), w.cend(), r2)));
    node u = G.randomNeighbor(r);
    node w;
    do {
      w = G.randomNeighbor(r);
    } while (w == u);
    if(G.hasEdge(u, w)) {
      l++;
    }
  }

  return static_cast<double>(l) / k;
}

} /* namespace NetworKit */
