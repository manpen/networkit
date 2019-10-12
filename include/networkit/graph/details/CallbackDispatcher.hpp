/**
 * CallbackDispatcher.hpp
 */

// networkit-format

#ifndef NETWORKIT_GRAPH_DETAILS_EDGE_LAMBDAS_HPP_
#define NETWORKIT_GRAPH_DETAILS_EDGE_LAMBDAS_HPP_

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>

namespace NetworKit {
namespace GraphDetails {

/*
 * In the following definition, Aux::FunctionTraits is used in order to only
 * execute lambda functions with the appropriate parameters. The
 * decltype-return type is used for determining the return type of the
 * lambda (needed for summation) but also determines if the lambda accepts
 * the correct number of parameters. Otherwise the return type declaration
 * fails and the function is excluded from overload resolution. Then there
 * are multiple possible lambdas with three (third parameter id or weight)
 * and two (second parameter can be second node id or edge weight for
 * neighbor iterators). This is checked using Aux::FunctionTraits and
 * std::enable_if. std::enable_if only defines the type member when the
 * given bool is true, this bool comes from std::is_same which compares two
 * types. The function traits give either the parameter type or if it is out
 * of bounds they define type as void.
 */

/**
 * Triggers a static assert error when no other method is chosen. Because of
 * the use of "..." as arguments, the priority of this method is lower than
 * the priority of the other methods. This method avoids ugly and unreadable
 * template substitution error messages from the other declarations.
 */
template <class F, void * = (void *)0>
typename Aux::FunctionTraits<F>::result_type callEdgeHandler(F &, ...) {
  // the strange condition is used in order to delay the evaluation of the
  // static assert to the moment when this function is actually used
  static_assert(!std::is_same<F, F>::value,
                "Your lambda does not support the required parameters or the "
                "parameters have the wrong type.");
  return std::declval<typename Aux::FunctionTraits<
      F>::result_type>(); // use the correct return type (this won't
                          // compile)
}

/**
 * Calls the given function f if its fourth argument is of the type edgeid
 * and third of type edgeweight Note that the decltype check is not enough
 * as edgeweight can be casted to node and we want to assure that .
 */
template <
    class F,
    typename std::enable_if<
        (Aux::FunctionTraits<F>::arity >= 3) &&
        std::is_same<edgeweight, typename Aux::FunctionTraits<F>::template arg<
                                     2>::type>::value &&
        std::is_same<edgeid, typename Aux::FunctionTraits<F>::template arg<
                                 3>::type>::value>::type * = (void *)0>
auto callEdgeHandler(F &f, node u, node v, edgeweight ew, edgeid id)
    -> decltype(f(u, v, ew, id)) {
  return f(u, v, ew, id);
}

/**
 * Calls the given function f if its third argument is of the type edgeid,
 * discards the edge weight Note that the decltype check is not enough as
 * edgeweight can be casted to node.
 */
template <class F,
          typename std::enable_if<
              (Aux::FunctionTraits<F>::arity >= 2) &&
              std::is_same<edgeid, typename Aux::FunctionTraits<
                                       F>::template arg<2>::type>::value &&
              std::is_same<node, typename Aux::FunctionTraits<F>::template arg<
                                     1>::type>::value /* prevent f(v, weight,
                                                       * eid)
                                                       */
              >::type * = (void *)0>
auto callEdgeHandler(F &f, node u, node v, edgeweight, edgeid id)
    -> decltype(f(u, v, id)) {
  return f(u, v, id);
}

/**
 * Calls the given function f if its third argument is of type edgeweight,
 * discards the edge id Note that the decltype check is not enough as node
 * can be casted to edgeweight.
 */
template <
    class F,
    typename std::enable_if<
        (Aux::FunctionTraits<F>::arity >= 2) &&
        std::is_same<edgeweight, typename Aux::FunctionTraits<F>::template arg<
                                     2>::type>::value>::type * = (void *)0>
auto callEdgeHandler(F &f, node u, node v, edgeweight ew, edgeid /*id*/)
    -> decltype(f(u, v, ew)) {
  return f(u, v, ew);
}

/**
 * Calls the given function f if it has only two arguments and the second
 * argument is of type node, discards edge weight and id Note that the
 * decltype check is not enough as edgeweight can be casted to node.
 */
template <class F,
          typename std::enable_if<
              (Aux::FunctionTraits<F>::arity >= 1) &&
              std::is_same<node, typename Aux::FunctionTraits<F>::template arg<
                                     1>::type>::value>::type * = (void *)0>
auto callEdgeHandler(F &f, node u, node v, edgeweight /*ew*/, edgeid /*id*/)
    -> decltype(f(u, v)) {
  return f(u, v);
}

/**
 * Calls the given function f if it has only two arguments and the second
 * argument is of type edgeweight, discards the first node and the edge id
 * Note that the decltype check is not enough as edgeweight can be casted to
 * node.
 */
template <
    class F,
    typename std::enable_if<
        (Aux::FunctionTraits<F>::arity >= 1) &&
        std::is_same<edgeweight, typename Aux::FunctionTraits<F>::template arg<
                                     1>::type>::value>::type * = (void *)0>
auto callEdgeHandler(F &f, node, node v, edgeweight ew, edgeid /*id*/)
    -> decltype(f(v, ew)) {
  return f(v, ew);
}

/**
 * Calls the given function f if it has only one argument, discards the
 * first node id, the edge weight and the edge id
 */
template <class F, void * = (void *)0>
auto callEdgeHandler(F &f, node, node v, edgeweight, edgeid) -> decltype(f(v)) {
  return f(v);
}

/**
 * Calls the given BFS handle with distance parameter
 */
template <class F>
auto callBFSHandler(F &f, node u, count dist) -> decltype(f(u, dist)) {
  return f(u, dist);
}

/**
 * Calls the given BFS handle without distance parameter
 */
template <class F> auto callBFSHandler(F &f, node u, count) -> decltype(f(u)) {
  return f(u);
}

} // namespace GraphDetails
} // namespace NetworKit

#endif // NETWORKIT_GRAPH_DETAILS_EDGE_LAMBDAS_HPP_
