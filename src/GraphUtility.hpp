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

#ifndef __PAP_BCO_PARSER__GRAPH_UTILITY__HPP
#define __PAP_BCO_PARSER__GRAPH_UTILITY__HPP

#include <queue>
#include <vector>
#include <boost/graph/two_bit_color_map.hpp>

namespace pap_solver {

class GraphUtility {
 public:
  template<typename Graph>
  static bool graph_contains_loop(const Graph& graph);

  /// @brief Prints the graph pass as an adjacency matrix.
  ///
  /// @param [in] graph        The matrix you want to print on the stream.
  /// @param [out] out_stream  The output stream where to print the adjancency
  ///                          matrix.
  /// @param [in] header       'True' if you want to print the header
  ///                          which display the num of vertices and edges.
  /// @param separator         The seprator character.
  ///
  template<typename Graph>
  static void printGraph_as_adjacencyMatrix(const Graph& graph,
                                            std::ostream* out_stream,
                                            bool header = false,
                                            char separator = ',');

  /// @brief A Vertex filter.
  typedef std::vector<bool> VertexFilter;

  /// @brief Discovers all disjointed graph (subgraph) in the graph passed
  ///        as paramter.
  /// @param [in] graph             The main graph which could contains
  ///                               disjointed sub graphs.
  /// @param [out] dis_graphs       A vector of vertex filter. Each vertex
  ///                               filter represents a map which says
  ///                               whether a vertex is present or not.
  /// @param [in] unitary_subgraph  If true even unitary node will be
  ///                               considered as sub graphs.
  /// @note It's recommended to avoid unitary_subgraph because of perfomance
  ///       reasons. Less memory will be used.
  /// @note The complexity should be O(|V|).
  template<typename Graph>
  static void find_all_disjointed_graph(
      const Graph& graph,
      std::vector<VertexFilter>* dis_graphs,
      bool unitary_subgraph = false);
};

template<typename Graph>
void GraphUtility::find_all_disjointed_graph(
    const Graph& graph,
    std::vector<VertexFilter>* dis_graphs,
    bool unitary_subgraph) {
  assert(dis_graphs != nullptr);

  // Some type definitions
  typedef typename Graph::vertex_descriptor VertexType;
  typedef std::queue<VertexType> OpenList;
  typedef std::vector<bool> ClosedList;

  // Some local variable
  const auto num_vertices = boost::num_vertices(graph);
  VertexFilter* current_filter = nullptr;

  // Clean the output
  dis_graphs->clear();

  // Initialize the closed list. All vertices are unexplored
  ClosedList closedList(num_vertices, false);
  size_t explored_node = 0;

  // Inizialize the open list. Insert the first vertex
  OpenList openList;
  const auto first_vertex = boost::vertices(graph).first;
  openList.push(*first_vertex);

  while (explored_node < num_vertices) {
    // Check if the open list is empty. In that case find another node
    if (openList.empty() == true) {
      // Force next iteration to create a new filter
      current_filter = nullptr;

      const auto all_vertices = boost::vertices(graph);
      const auto finder = std::find_if(
          all_vertices.first, all_vertices.second,
          [&] (const VertexType& v) {
            return (closedList[v] == false);
          });

      // There must to be another vertex to explore
      assert(finder != all_vertices.second);

      // Insert this new node in the openlist
      openList.push(*finder);
    }

    while (openList.empty() == false) {
      // Get the first vertex in the open list and remove from it
      const auto current_vertex = openList.front();
      openList.pop();

      // Check if the current vertex is a isolated node.
      // In that case, and the option for unitary graph is disabled
      // then you can skip that node
      if (unitary_subgraph == true ||
          boost::out_degree(current_vertex, graph) > 0) {
        // If the current filter is null, create a new one and
        // use it
        if (current_filter == nullptr) {
          dis_graphs->emplace_back(num_vertices, false);
          current_filter = &(dis_graphs->back());
        }

        // Add the current vertex to the current filter
        (*current_filter)[current_vertex] = true;

        // Add it to the closed list
        closedList[current_vertex] = true;
        ++explored_node;

        // Get all adjacent vertices from the current vertex
        const auto adj_vertices = boost::adjacent_vertices(current_vertex,
                                                           graph);

        // Add all adjacent vertices to the open list
        std::for_each(
            adj_vertices.first, adj_vertices.second,
            [&] (const VertexType& a_v) {
              // Check the vertex is not in the closed list
              if (closedList[a_v] == false) {
                openList.push(a_v);
              }
            });
      }
    }  // end while on openList
  }  // end while on all vertices explored
}

template<typename Graph>
void GraphUtility::printGraph_as_adjacencyMatrix(const Graph& graph,
                                                 std::ostream* out_stream,
                                                 bool header,
                                                 char separator) {
  assert(out_stream != nullptr);

  // Some type definitions
  typedef typename Graph::vertex_descriptor VertexType;
  typedef typename Graph::edge_descriptor EdgeType;

  // Const definition and local variables
  const auto num_vertices = boost::num_vertices(graph);
  const auto num_edges = boost::num_edges(graph);
  auto& os = *out_stream;

  // Print header if specified in parameters
  if (header) {
    os << num_vertices << ',' << num_edges << '\n';
  }

  // Loop on all vertices O(|V|)
  for (size_t i = 0; i < num_vertices; ++i) {
    auto adj_vertices = boost::adjacent_vertices(i, graph);
    std::vector<char> row_adj_vertex(num_vertices, '0');

    // Loop on all adjacent vertices of i-th vertex
    std::for_each(
        adj_vertices.first, adj_vertices.second,
        [&] (const VertexType& adj_v) {
          row_adj_vertex[adj_v] = '1';
        });

    // Print the row vector on the output stream
    size_t col = 0;
    for (const auto& c : row_adj_vertex) {
      os << c << (++col < num_vertices ? separator : '\n');
    }
  }
}

template<typename Graph>
bool GraphUtility::graph_contains_loop(const Graph& graph) {
  // Type definitions
  typedef typename Graph::vertex_descriptor VertexType;
  typedef typename Graph::edge_descriptor EdgeType;
  typedef boost::two_bit_color_map<> color_map_t;
  typedef std::queue<VertexType> open_list_t;

  // Color definitions
  static const auto white_t =  boost::color_traits<
    boost::two_bit_color_type>::white();
  static const auto black_t =  boost::color_traits<
    boost::two_bit_color_type>::black();
  static const auto gray_t =  boost::color_traits<
    boost::two_bit_color_type>::gray();

  const auto num_vertices = boost::num_vertices(graph);
  color_map_t color_map(num_vertices);

  open_list_t openlist;
  const auto root = *(vertices(graph).first);
  openlist.push(root);
  put(color_map, root, black_t);

  while (openlist.empty() == false) {
    const auto& node = openlist.front();
    auto this_node_color = get(color_map, node);
    auto edges_list = out_edges(node, graph);
    for (auto i = edges_list.first;
         i != edges_list.second;
         ++i) {
      auto target = boost::target(*i, graph);
      auto color_adjacent = get(color_map, target);
      if (color_adjacent == white_t) {
        if (this_node_color == black_t)
          put(color_map, target, gray_t);
        else
          put(color_map, target, black_t);
        openlist.push(target);
      } else if (color_adjacent == this_node_color) {
        return true;
      }
    }
    openlist.pop();
  }
  return false;
}

}  // namespace pap_solver

#endif  // __PAP_BCO_PARSER__GRAPH_UTILITY__HPP
