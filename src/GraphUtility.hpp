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
#include <stack>
#include <vector>
#include <utility>
#include <istream>
#include <ostream>
#include <boost/graph/random.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
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

  static bool check_all_disjointed(
      const std::vector<VertexFilter>& dis_graphs,
      const size_t total_num_vertices) noexcept;

  /// @brief Generate a random graph. Usefull for test purpose.
  ///
  /// @param [in] num_vertices   The number of vertices you want in the graph.
  /// @param [in] num_edges      The number of edges you want in the graph.
  /// @param [out] output_graph  The output graph.
  /// @param [in,out] rnd_engine A random engine.
  ///
  /// @note Actually the number of of the output graph edges may be different
  /// that what specified in the parameter. That because additional edges
  /// will be added in order to fix a disconnected graph.
  template<typename Graph, typename RndEngine>
  static void generate_random_graph(const size_t num_vertices,
                                    const size_t num_edges,
                                    Graph* output_graph,
                                    RndEngine* rnd_engine);

  template<typename Graph>
  static void printGraph_asArchive(const Graph& graph,
                                   std::ostream* out_stream);

  template<typename Graph>
  static void readGraph_asArchive(std::istream* in_stream, Graph* graph);
};

template<typename Graph>
void GraphUtility::printGraph_asArchive(const Graph& graph,
                                        std::ostream* out_stream) {
  assert(out_stream != nullptr);

  boost::archive::text_oarchive archive(*out_stream);
  archive << graph;
}

template<typename Graph>
void GraphUtility::readGraph_asArchive(std::istream* in_stream,
                                       Graph* graph) {
  assert(in_stream != nullptr);
  assert(graph != nullptr);

  boost::archive::text_iarchive archive(*in_stream);
  archive >> *graph;
}

template<typename Graph, typename RndEngine>
void GraphUtility::generate_random_graph(const size_t num_vertices,
                                         const size_t num_edges,
                                         Graph* output_graph,
                                         RndEngine* rnd_engine) {
  assert(output_graph != nullptr);
  assert(rnd_engine != nullptr);

  static_assert(std::is_integral<typename Graph::vertex_descriptor>::value,
                "The vertex type (descriptor) must to be a integer type!");

  // Some local variables
  std::vector<VertexFilter> sub_graphs;

  // Clean the output
  output_graph->clear();

  // Generate the random graph
  boost::generate_random_graph(*output_graph, num_vertices, num_edges,
                               *rnd_engine, false);

  // The graph could be disjointed, check that
  find_all_disjointed_graph(*output_graph, &sub_graphs, true);

  // While there are disjointed graph
  size_t num_sub_graphs;
  while ((num_sub_graphs = sub_graphs.size()) > 1) {
    // Pick two of those (rndly) graph and join them
    std::uniform_int_distribution<int> rnd_sub_graph(0, num_sub_graphs - 1);

    // TODO(biagio): no efficient solution to generate different numbers
    auto index1 = rnd_sub_graph(*rnd_engine);
    decltype(index1) index2;
    while ((index2 = rnd_sub_graph(*rnd_engine)) == index1) { }

    // Get those two graph
    const auto& graph1 = sub_graphs[index1];
    const auto& graph2 = sub_graphs[index2];

    // Pick two different vertices
    int vertex1 = -1;
    int vertex2 = -1;
    for (size_t i = 0;
         i < num_vertices && (vertex1 == -1 || vertex2 == -1);
         ++i) {
      if (vertex1 == -1 && graph1[i] == true) {
        vertex1 = i;
      }
      if (vertex2 == -1 &&  graph2[i] == true) {
        vertex2 = i;
      }
    }  // for
    assert(vertex1 != vertex2);
    assert(graph1[vertex1] == true);
    assert(graph2[vertex2] == true);

    // Connect those two vertices
    boost::add_edge(vertex1, vertex2, *output_graph);

    // Erase one of those graph
    sub_graphs.erase(sub_graphs.cbegin() + index2);

    // TODO(biagio): l'algoritmo non è proprio ottimale.
    // Innanzittutto i vertici sono selezionati come i primi che trova mentre
    // potrebbero essere pickati randomicamente
    // Inoltre dovresti aggiornare il subgrafo1 mettendo a true
    // i nuovi vertici importati ora da subgraph2.
    // Non mettendoli, come in questo caso, alla prossima iterazione
    // i vertici di graph2 (ora in graph1) non vengono visti e quindi
    // non selezionabili.
  }  // end while
}

bool GraphUtility::check_all_disjointed(
    const std::vector<VertexFilter>& dis_graphs,
    const size_t total_num_vertices) noexcept {
  std::vector<int> vertex_in_graphID(total_num_vertices, -1);

  size_t current_index_subgraph = 0;

  // Cycle on all sub graphs
  for (const auto& g : dis_graphs) {
    size_t current_index_vertex = 0;

    // Cycle on all vertices of this sub graph
    for (const auto& v : g) {
      if (v == true) {
        auto& vertex_belogs = vertex_in_graphID[current_index_vertex];
        if (vertex_belogs != -1) {
          return false;
        } else {
          vertex_belogs = current_index_subgraph;
        }
      }

      // Incremente the index of vertex
      ++current_index_vertex;
    }

    // Increment the index of sub graph
    ++current_index_subgraph;
  }

  return true;
}

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
  ++explored_node;
  closedList[*first_vertex] = true;

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
      closedList[*finder] = true;
      ++explored_node;
    }

    while (openList.empty() == false) {
      // Get the first vertex in the open list and remove from it
      const auto current_vertex = openList.front();
      openList.pop();

      // Check if the current vertex is a isolated node.
      // In that case, and the option for unitary graph is disabled
      // then you can skip that node
      if (unitary_subgraph == true ||
          (boost::out_degree(current_vertex, graph) > 0)) {
        // If the current filter is null, create a new one and
        // use it
        if (current_filter == nullptr) {
          dis_graphs->emplace_back(num_vertices, false);
          current_filter = &(dis_graphs->back());
        }

        // Add the current vertex to the current filter
        (*current_filter)[current_vertex] = true;

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

                // Add it to the closed list
                closedList[a_v] = true;
                ++explored_node;
              }
            });
      }
    }  // end while on openList
  }  // end while on all vertices explored

  // Assertion on disjointed property
  assert(GraphUtility::check_all_disjointed(
      *dis_graphs, num_vertices) == true);
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
  typedef std::stack<std::pair<VertexType, VertexType>> OpenList;
  typedef std::vector<bool> ClosedList;

  // Some local variables
  size_t num_node_explored = 0;
  const auto num_vertices = boost::num_vertices(graph);
  if (num_vertices == 0) return false;

  // Closed List inizialization
  ClosedList closedList(num_vertices, false);

  // Open List inizialization
  OpenList openList;
  const auto p_firstVertex = boost::vertices(graph).first;
  const auto p_nullVertex = boost::vertices(graph).second;
  openList.push(std::make_pair(*p_firstVertex, *p_nullVertex));
  closedList[*p_firstVertex] = true;
  ++num_node_explored;

  while (openList.empty() == false) {
    // Pop the next node
    const auto current_element = openList.top();
    const auto& current_node = current_element.first;
    const auto& current_parent = current_element.second;
    openList.pop();

    // Explore all adjacent vertices
    const auto adj_vertices = boost::adjacent_vertices(current_node, graph);

    // Look for a back vertices which is not the parent
    const auto finder = std::find_if(
        adj_vertices.first, adj_vertices.second,
        [&] (const VertexType& a_v) {
          if (closedList[a_v] == true) {
            if (a_v != current_parent) {
              // Found a back edge. A cycle is here!
              return true;
            }
          } else {
            openList.push(std::make_pair(a_v, current_node));
            closedList[a_v] = true;
            ++num_node_explored;
          }
          return false;
        });

    if (finder != adj_vertices.second) {
      return true;
    }
  }  // while openList is not empty

  // Assertion all nodes are explored. In case the graph is disconnected
  assert(num_node_explored == num_vertices);

  // No cycle found!
  return false;
}

}  // namespace pap_solver

#endif  // __PAP_BCO_PARSER__GRAPH_UTILITY__HPP
