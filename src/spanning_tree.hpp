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

  /// An associative container for each vertex it returns its parent.
  typedef std::map<VertexType, VertexType> MapVertexParent;

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


  /// A null vertex. The parent of the root will be equal to this.
  const VertexType sNullVertex = Graph::null_vertex();

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
  void print(std::ostream* os) const;

  /// @brief Clean the tree.
  void clear() noexcept;

  inline const MapVertexParent& get_mapped_spanning_tree() const noexcept;

  inline const EdgeFilter& get_edges_in_spanning_tree() const noexcept;

  inline const FilteredGraph get_filtered_graph(const Graph& graph);

 private:
  MapVertexParent m_mapped_spanning_tree;
  EdgeFilter m_edges_filter;

  template<typename RND>
  static void makeMappedFromGraph(const Graph& graph,
                                  RND* rnd_engine,
                                  MapVertexParent* output_map);

  static void makeFilterFromMap(const Graph& graph,
                                const MapVertexParent& input_map,
                                EdgeFilter* output_filter);
};

template <typename Graph>
template<typename RandomEngine>
void SpanningTree<Graph>::generate_rnd_spanning_tree(const Graph& graph,
                                                     RandomEngine* rnd_engine) {
  assert(rnd_engine != nullptr);

  clear();
  makeMappedFromGraph(graph, rnd_engine, &m_mapped_spanning_tree);
  makeFilterFromMap(graph, m_mapped_spanning_tree, &m_edges_filter);
}

template<typename Graph>
void SpanningTree<Graph>::clear() noexcept {
  m_mapped_spanning_tree.clear();
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
void SpanningTree<Graph>::print(std::ostream* os) const {
  assert(os != nullptr);
  // TODO(biagio): to implement
}

template<typename Graph>
inline
const typename SpanningTree<Graph>::MapVertexParent&
SpanningTree<Graph>::get_mapped_spanning_tree() const noexcept {
  return m_mapped_spanning_tree;
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


}  // namespace pap_solver

#endif

