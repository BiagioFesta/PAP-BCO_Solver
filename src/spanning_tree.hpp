// Copyright 2016 <Biagio Festa>

/*
  This file is part of PAP-BCO_Solver.

  PAP-BCO_Solver is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
    
  PAP-BCO_Solver is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
    
  You should have received a copy of the GNU General Public License
  along with PAP-BCO_Solver.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __PAP_BCO_SOLVER__SPANNING_TREE__HPP
#define __PAP_BCO_SOLVER__SPANNING_TREE__HPP

#include <map>
#include <ostream>
#include <type_traits>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>
#include "GraphUtility.hpp"

namespace pap_solver {

/// @brief Spanning Tree for a Graph
/// @author Biagio Festa
/// @template Graph must to be a concept of boost Graph
template<typename Graph>
class SpanningTree {
 public:
  /// A vertex type of the graph.
  typedef typename Graph::vertex_descriptor VertexType;
  /// A edge type of the graph.
  typedef typename Graph::edge_descriptor EdgeType;

  /// An associative container for each edges it returns if in spanning tree.
  typedef std::map<EdgeType, bool> EdgeFilter;

  /// @brief A predicate (function object) which return
  ///        whether an edge belongs to the filtered graph or not.
  struct PredicateFilterEdge {
    PredicateFilterEdge() = default;

    explicit PredicateFilterEdge(EdgeFilter* map) :
        m_edges_map(*map) {
    }
    inline bool operator()(const EdgeType e) const {
      return boost::get(m_edges_map, e);
    }

    boost::associative_property_map<EdgeFilter> m_edges_map;
  };

  /// A view on a graph.
  typedef boost::filtered_graph<Graph, PredicateFilterEdge> FilteredGraph;

  static_assert(std::is_integral<VertexType>::value,
                "The vertex type (descriptor) must to be a integer type!");

  /// @brief Default constructor.
  SpanningTree() = default;

  /// @brief Default destructor.
  ~SpanningTree() = default;

  template<typename RandomEngine>
  void generate_rnd_spanning_tree(const Graph& graph,
                                  RandomEngine* rnd_engine);

  /// @brief Print the subtree with the notation:
  ///                   (VERTEX)  ->  (vertex PARENT)
  ///                             ...
  ///
  /// @param os [out]   The output stream wehere the subtree
  ///                   will be printed.
  ///
  void print(const Graph& graph, std::ostream* os) const;

  /// @brief Clean the tree.
  void clear() noexcept;

  void perform_transformation(const Graph& graph,
                              const EdgeType& edge_to_add,
                              const EdgeType& edge_to_remove);

  inline const EdgeFilter& get_edges_in_spanning_tree() const noexcept;

  inline const FilteredGraph get_filtered_graph(const Graph& graph);

 private:
  typedef std::map<VertexType, VertexType> MapVertexParent;

  EdgeFilter m_edges_filter;

  template<typename RND>
  static void makeMappedFromGraph(const Graph& graph,
                                  RND* rnd_engine,
                                  MapVertexParent* output_map);

  static void makeFilterFromMap(const Graph& graph,
                                const MapVertexParent& input_map,
                                EdgeFilter* output_filter);

  bool this_is_a_valid_spanning_tree(const Graph& graph);
};

template <typename Graph>
template<typename RandomEngine>
void SpanningTree<Graph>::generate_rnd_spanning_tree(const Graph& graph,
                                                     RandomEngine* rnd_engine) {
  assert(rnd_engine != nullptr);

  clear();
  MapVertexParent parent_map;
  makeMappedFromGraph(graph, rnd_engine, &parent_map);
  makeFilterFromMap(graph, parent_map, &m_edges_filter);

  assert(this_is_a_valid_spanning_tree(graph) == true);
}

template<typename Graph>
void SpanningTree<Graph>::clear() noexcept {
  m_edges_filter.clear();
}

template<typename Graph>
template<typename RND>
void SpanningTree<Graph>::makeMappedFromGraph(const Graph& graph,
                                              RND* rnd_engine,
                                              MapVertexParent* output_map) {
  assert(output_map != nullptr && rnd_engine != nullptr);

  output_map->clear();
  if (boost::num_vertices(graph) == 0) {
    throw std::runtime_error("The graph has no vertices!");
  }

  boost::associative_property_map<MapVertexParent> predecessor_map(*output_map);
  boost::random_spanning_tree(graph,
                              *rnd_engine,
                              boost::predecessor_map(predecessor_map).
                              vertex_index_map(boost::identity_property_map()));
}

template<typename Graph>
void SpanningTree<Graph>::makeFilterFromMap(const Graph& graph,
                                            const MapVertexParent& input_map,
                                            EdgeFilter* output_filter) {
  assert(output_filter != nullptr);

  auto& edges_in_spanning_tree = *output_filter;

  // Initialization (all false)
  output_filter->clear();
  auto edges_in_graph = edges(graph);
  std::for_each(edges_in_graph.first,
                edges_in_graph.second,
                [&edges_in_spanning_tree]
                (const EdgeType& e) {
                  edges_in_spanning_tree[e] = false;
                });

  // TODO(biagio): could be improved
  for (const auto& node : input_map) {
    const auto edges_from_this = out_edges(node.first, graph);
    std::for_each(edges_from_this.first,
                  edges_from_this.second,
                  [&edges_in_spanning_tree, &node, &graph]
                  (const EdgeType& e) {
                    if (target(e, graph) == node.second)
                      edges_in_spanning_tree[e] = true;
                  });
  }
}

template<typename Graph>
void SpanningTree<Graph>::print(const Graph& graph, std::ostream* os) const {
  if (os == nullptr) return;

  const auto verxs = vertices(graph);
  std::for_each(verxs.first,
                verxs.second,
                [this, &os, &graph]
                (const VertexType& v) {
                  if (v == 0) {
                    *os << '(' << v << ") ---> (ROOT)\n";
                  } else {
                    const auto& links = out_edges(v, graph);
                    auto finder = std::find_if(links.first,
                                               links.second,
                                               [this]
                                               (const EdgeType& e) {
                                                 return m_edges_filter.at(e);
                                               });

                    if (finder != links.second) {
                      *os << '(' << v << ") ---> (" <<
                          target(*finder, graph) << ")\n";
                    }
                  }
                });
}

template<typename Graph>
inline
const typename SpanningTree<Graph>::EdgeFilter&
SpanningTree<Graph>::get_edges_in_spanning_tree() const noexcept {
  return m_edges_filter;
}

template<typename Graph>
inline
const typename SpanningTree<Graph>::FilteredGraph
SpanningTree<Graph>::get_filtered_graph(const Graph& graph) {
  return FilteredGraph(graph, PredicateFilterEdge(&m_edges_filter));
}

template<typename Graph>
bool SpanningTree<Graph>::this_is_a_valid_spanning_tree(
    const Graph& graph) {
  // A valid spanning tree must have no cycle
  const auto filtered_graph = get_filtered_graph(graph);

  // Check whether the spanning tree contains loops.
  return GraphUtility::graph_contains_loop(filtered_graph) ?
      false : true;
}

template<typename Graph>
void SpanningTree<Graph>::perform_transformation(
    const Graph& graph, const EdgeType& edge_to_add,
    const EdgeType& edge_to_remove) {
#ifdef _DEBUG
  // Check condition od edge
  assert(m_edges_filter.at(edge_to_remove) == true);
  assert(m_edges_filter.at(edge_to_add) == false);
#endif

  m_edges_filter[edge_to_remove] = false;
  m_edges_filter[edge_to_add] = true;
}


}  // namespace pap_solver

#endif

