/*
 * UDegreeIndex.cpp
 *
 *  Created on: 01.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <networkit/linkprediction/UDegreeIndex.hpp>

namespace NetworKit {

double UDegreeIndex::runImpl(node u, node) {
  return static_cast<double>(G->degree(u));
}

} // namespace NetworKit
