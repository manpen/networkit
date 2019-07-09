/*
 * DotGraphWriter.hpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#ifndef DOTGRAPHWRITER_H
#define DOTGRAPHWRITER_H

#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphWriter.hpp>

namespace NetworKit {

/**

 * @ingroup io
 *
 * This class turns a graph into a very basic GraphViz file as documented in the official manual [1].
 * If a more thorough support is desired, please contact the developers: https://github.com/networkit/networkit
 *
 * [1] https://graphviz.gitlab.io/_pages/pdf/dotguide.pdf
 */
class DotGraphWriter final : public GraphWriter {
public:
	/**
	 * Write a graph as a GraphViz/file.
	 *
	 * @param[in]	graph	The graph object
	 * @param[in]	path	The file path to be written to
	 */
	void write(const Graph &G, const std::string &path) const override;

};

} /* namespace NetworKit */
#endif /* DOTGRAPHWRITER_H */
