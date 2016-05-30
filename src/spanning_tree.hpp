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

  static_assert(std::is_integral<VertexType>::value,
                "The vertex type (descriptor) must to be a integer type!");

  /// @brief Default constructor.
  SpanningTree() noexcept;

  /// @brief Create a spanning tree with a random approach.
  ///
  /// @param [in] g           The starting graph.
  /// @param rnd_engine       A valid random engine.
  ///
  /// @note The graph must to have all edges with unitary cost.
  /// @note The root of the tree will be the first vertex of the graph.
  template<typename RND>
  void makeRandom_fromGraph(const Graph& g, RND* rnd_engine);

  /// @brief Print the subtree with the notation:
  ///                   (VERTEX)  ->  (vertex PARENT)
  ///                             ...
  ///
  /// @param os [out]   The output stream wehere the subtree
  ///                   will be printed.
  ///
  void print(std::ostream* os) const;

 private:
  /// An associative container with <key,value> VertexType.
  typedef std::map<VertexType, VertexType> Map;

  Map m_data_tree;

  /// A null vertex. The parent of the root will be equal to this.
  VertexType sNullVertex;
};

template<typename Graph>
SpanningTree<Graph>::SpanningTree() noexcept:
                                     sNullVertex(Graph::null_vertex()) {
}

template<typename Graph>
template<typename RND>
void SpanningTree<Graph>::makeRandom_fromGraph(const Graph& g,
                                               RND* rnd_engine) {
  boost::associative_property_map<Map> predecessor_map(m_data_tree);
  boost::random_spanning_tree(g,
                              *rnd_engine,
                              boost::predecessor_map(predecessor_map).
                              vertex_index_map(boost::identity_property_map()));
  if (m_data_tree.begin()->second != sNullVertex) {
    throw std::runtime_error("I cannot create a proper spanning tree. "
                             "The first vertex is not the root!");
  }
}

template<typename Graph>
void SpanningTree<Graph>::print(std::ostream* os) const {
  if (os == nullptr || m_data_tree.size() == 0) return;

  for (const auto& node : m_data_tree) {
    if (node.second == sNullVertex) {
      *os << node. first << " -> " << "[ROOT]" << '\n';
    } else {
      *os << node.first << "  ->  " << node.second << '\n';
    }
  }
}

}  // namespace pap_solver

#endif

