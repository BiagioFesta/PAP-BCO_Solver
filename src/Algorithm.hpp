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

#ifndef __PAP_BCO_PARSER__ALGORITHM__HPP
#define __PAP_BCO_PARSER__ALGORITHM__HPP

#include <map>
#include <algorithm>
#include <queue>
#include <ostream>
#include <utility>
#include <chrono>
#include "spanning_tree.hpp"
#include <boost/graph/filtered_graph.hpp>


namespace pap_solver {

template<typename Graph, typename RndGenerator>
class Algorithm {
 public:
  /// @brief A port assignment for each vertex of the graph
  enum PortAssignment {
    UnAssigned,
    PortA,
    PortB,
    PortAB
  };

  /// @brief Vertex Descriptor
  typedef typename Graph::vertex_descriptor VertexType;

  /// @brief Edge Descriptor
  typedef typename Graph::edge_descriptor EdgeType;

  /// @brief A predicate (function object) which return
  ///        whether an edge belongs to the filtered graph or not.
  template<typename AssociativeMap>
  struct PredicateFilterEdge {
    typedef boost::associative_property_map<AssociativeMap> edge_list_t;

    PredicateFilterEdge() = default;
    explicit PredicateFilterEdge(AssociativeMap* map) :
        m_edges_list(*map) {
    }
    bool operator()(const EdgeType e) const {
      return boost::get(m_edges_list, e);
    }
    edge_list_t m_edges_list;
  };

  /// @brief An map assignment.
  typedef std::map<VertexType, PortAssignment> MapAssignment;

  /// @brief A filter on edges. Whether they belong to the filtered
  ///        graph or no.
  typedef std::map<EdgeType, bool> EdgesMapFilter;

  /// @brief A struct with the solution and some its details.
  struct Solution {
    MapAssignment m_assignment;
    size_t m_size_solution;
    SpanningTree<Graph> m_mapped_spanning_tree;
    EdgesMapFilter m_edges_into_spanning_tree;
    std::chrono::milliseconds m_time_for_solution;
  };

  /// @brief Default constructor.
  Algorithm() = default;

  bool solve(const Graph& graph, Solution* out_solution);

  /// @brief Sets the seed for the random engine.
  void set_seed(int seed);

  static void print_solution(std::ostream* os,
                             const Solution& solution);

  static void print_spanning_tree(std::ostream* os,
                                  const Solution& solution);


 private:
  RndGenerator m_rnd_engine;

  /// @brief Solve the port assignment problem.
  ///
  /// @param [in] graph          A valid graph.
  /// @param [in] spanning_tree  A valid spanning tree from that graph.
  /// @param [in] tree_map       A filter which describes which vertices
  ///                            belong to the spanning tree.
  /// @param [out] out_solution  It fills the solution's fields assignment and
  ///                            size of solution.
  template<typename SpanningTreeT>
  void solve_problem(const Graph& graph,
                     const SpanningTreeT& spanning_tree,
                     const EdgesMapFilter& tree_map,
                     Solution* out_solution) const;

  /// @param [in] graph       A valid graph.
  /// @param [in] tree_map    A filter which described which vertices
  ///                         belong to the spanning tree.
  /// @param [in] e_to_test   The edge you want to test.
  ///
  /// @return 'true' if the edge is a odd co-tree edge, 'false'
  //          otherwise.
  /// @note 'e_to_test' MUST to be an co-tree edge. The check won't
  ///        performed by the function!
  bool is_odd_cotree_edge(const Graph& graph,
                          const EdgesMapFilter& tree_map,
                          const EdgeType& e_to_test) const;

  /// @brief Generate a filter for a random spanning tree.
  /// @param [in] graph            A valid graph.
  /// @param [out] out_solution    It fills the solution's fields
  ///                              the mapped spanning tree and the filter.
  void generate_random_filter_tree(const Graph& graph,
                                   Solution* out_solution);

  template<typename SpanningTreeT>
  void assign_accordance_spanning_tree(const SpanningTreeT& spanning_tree,
                                       Solution* out_solution) const;
};

template<typename Graph, typename RndGenerator>
bool Algorithm<Graph, RndGenerator>::solve(
    const Graph& graph, Solution* out_solution) {
  typedef boost::filtered_graph<
    Graph, PredicateFilterEdge<EdgesMapFilter>> FilteredGraph;

  auto time_start = std::chrono::steady_clock::now();

  if (out_solution == nullptr) {
    auto time_stop = std::chrono::steady_clock::now();
    out_solution->m_time_for_solution =
        std::chrono::duration_cast<std::chrono::milliseconds>(
        time_stop - time_start);
    return false;
  }

  generate_random_filter_tree(graph, out_solution);
  EdgesMapFilter& spanning_tree_map =
      out_solution->m_edges_into_spanning_tree;

  PredicateFilterEdge<EdgesMapFilter>
      predicate_spanning_tree(&spanning_tree_map);
  FilteredGraph spanning_tree(graph, predicate_spanning_tree);

  solve_problem(graph, spanning_tree, spanning_tree_map, out_solution);

  auto time_stop = std::chrono::steady_clock::now();
  out_solution->m_time_for_solution =
      std::chrono::duration_cast<std::chrono::milliseconds>(
      time_stop - time_start);

  return true;
}

template<typename Graph, typename RndGenerator>
template<typename SpanningTreeT>
void Algorithm<Graph, RndGenerator>::solve_problem(
    const Graph& graph,
    const SpanningTreeT& spanning_tree,
    const EdgesMapFilter& tree_map,
    Solution* out_solution) const {
  using boost::vertices;
  using boost::edges;
  using boost::adjacent_vertices;

  // The value to return
  auto& rts = out_solution->m_assignment;

  // Preliminary assignment in accordance with the spanning tree
  assign_accordance_spanning_tree(spanning_tree, out_solution);

  // Find odd cotree edges
  EdgesMapFilter odd_cotree_edges;
  std::for_each(edges(graph).first,
                edges(graph).second,
                [this, &spanning_tree, &odd_cotree_edges, &graph, &tree_map]
                (const EdgeType& e) {
                  auto finder = tree_map.find(e);
                  if (finder == tree_map.cend() || finder->second == false) {
                    // 'e' is an co-tree edge, we're going to see
                    // whether it's odd or not.
                    if (is_odd_cotree_edge(graph, tree_map, e) == true) {
                      odd_cotree_edges[e] = true;
                    }
                  }
                });

  // Create a g' considering only odd cotree edges.
  PredicateFilterEdge<decltype(odd_cotree_edges)>
      predicate_odd_cotree_edges(&odd_cotree_edges);
  boost::filtered_graph<Graph, decltype(predicate_odd_cotree_edges)>
      g_prime(graph, predicate_odd_cotree_edges);

  bool found_one_degree;
  bool g_prime_empty = false;
  size_t number_of_AB = 0;
  while (!g_prime_empty) {
    found_one_degree = false;
    auto range_verts = vertices(g_prime);
    for (auto i = range_verts.first;
         i != range_verts.second && found_one_degree == false;
         ++i) {
      auto degree = out_degree(*i, g_prime);
      if (degree == 1) {
        EdgeType edge = *(out_edges(*i, g_prime).first);
        auto target = boost::target(edge, g_prime);
        rts[target] = PortAssignment::PortAB;
        ++number_of_AB;
        std::for_each(out_edges(target, g_prime).first,
                      out_edges(target, g_prime).second,
                      [&odd_cotree_edges]
                      (const EdgeType& e) {
                        odd_cotree_edges[e] = false;
                      });
        found_one_degree = true;
      }
    }
    if (found_one_degree == false) {
      size_t num_degree_max = 0;
      VertexType max_vertex_degree;
      for (auto i = vertices(g_prime).first;
           i != vertices(g_prime).second;
           ++i) {
        auto degree = out_degree(*i, g_prime);
        if (degree > num_degree_max) {
          num_degree_max = degree;
          max_vertex_degree = *i;
        }
      }
      if (num_degree_max > 0) {
        rts[max_vertex_degree] = PortAssignment::PortAB;
        ++number_of_AB;
        std::for_each(out_edges(max_vertex_degree, g_prime).first,
                      out_edges(max_vertex_degree, g_prime).second,
                      [&odd_cotree_edges]
                      (const EdgeType& e) {
                        odd_cotree_edges[e] = false;
                      });
      } else {
        g_prime_empty = true;
      }
    }
  }
  out_solution->m_size_solution = number_of_AB;
}

template<typename Graph, typename RndGenerator>
bool Algorithm<Graph, RndGenerator>::is_odd_cotree_edge(
    const Graph& graph,
    const EdgesMapFilter& tree_map,
    const EdgeType& e_to_test) const {
  using boost::vertices;
  using boost::edges;
  using boost::out_edges;
  using boost::get;
  using boost::put;

  EdgesMapFilter edges_in_spanning_tree_plus_e = tree_map;
  edges_in_spanning_tree_plus_e[e_to_test] = true;

  PredicateFilterEdge<decltype(edges_in_spanning_tree_plus_e)>
      predicate_spanning_tree_plus_e(&edges_in_spanning_tree_plus_e);
  boost::filtered_graph<Graph, decltype(predicate_spanning_tree_plus_e)>
      cycled_graph(graph, predicate_spanning_tree_plus_e);


  typedef boost::two_bit_color_map<> color_map_t;
  typedef std::queue<VertexType> open_list_t;
  static const auto white_t =  boost::color_traits<
    boost::two_bit_color_type>::white();
  static const auto black_t =  boost::color_traits<
    boost::two_bit_color_type>::black();
  static const auto gray_t =  boost::color_traits<
    boost::two_bit_color_type>::gray();

  const auto num_vertices = boost::num_vertices(cycled_graph);
  color_map_t color_map(num_vertices);

  open_list_t openlist;
  const auto root = *(vertices(cycled_graph).first);
  openlist.push(root);
  put(color_map, root, black_t);

  while (openlist.empty() == false) {
    const auto& node = openlist.front();
    auto this_node_color = get(color_map, node);
    auto edges_list = out_edges(node, cycled_graph);
    for (auto i = edges_list.first;
         i != edges_list.second;
         ++i) {
      auto target = boost::target(*i, cycled_graph);
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

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::generate_random_filter_tree(
    const Graph& graph, Solution* out_solution) {
  using boost::out_edges;
  using boost::target;
  using boost::edges;

  auto& mapped_spanning_tree = out_solution->m_mapped_spanning_tree;
  auto& edges_in_spanning_tree = out_solution->m_edges_into_spanning_tree;

  mapped_spanning_tree.clear();
  mapped_spanning_tree.makeRandom_fromGraph(graph, &m_rnd_engine);

  // Initialization of map property.
  edges_in_spanning_tree.clear();
  auto edges_in_graph = edges(graph);
  std::for_each(edges_in_graph.first,
                edges_in_graph.second,
                [&edges_in_spanning_tree]
                (const EdgeType& e) {
                  edges_in_spanning_tree[e] = false;
                });

  for (const auto& node : mapped_spanning_tree.getMap()) {
    std::for_each(out_edges(node.first, graph).first,
                  out_edges(node.first, graph).second,
                  [&edges_in_spanning_tree, &node, &graph]
                  (const EdgeType& e) {
                    if (target(e, graph) == node.second)
                      edges_in_spanning_tree[e] = true;
                  });
  }
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::set_seed(int seed) {
  m_rnd_engine.seed(seed);
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::print_solution(
    std::ostream* os,
    const Solution& solution) {
  if (os == nullptr) return;
  for (const auto& assignment : solution.m_assignment) {
    *os << "Vertex (" << assignment.first << ") ---> ";
    switch (assignment.second) {
      case PortAssignment::PortA:
        *os << "Port A";
        break;
      case PortAssignment::PortB:
        *os << "Port B";
        break;
      case PortAssignment::PortAB:
        *os << "Port AB";
        break;
      default:
        *os << "Unassigned!";
        break;
    }
    *os << '\n';
  }
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::print_spanning_tree(
    std::ostream* os,
    const Solution& solution) {
  solution.m_mapped_spanning_tree.print(os);
}

template<typename Graph, typename RndGenerator>
template<typename SpanningTreeT>
void Algorithm<Graph, RndGenerator>::assign_accordance_spanning_tree(
    const SpanningTreeT& spanning_tree,
    Solution* out_solution) const {
  using boost::vertices;
  using boost::adjacent_vertices;

  typedef std::queue<VertexType> open_list_t;

  auto& assignment = out_solution->m_assignment;

  // Initialization all vertices are unassigned O(V)
  const auto vertices_graph = vertices(spanning_tree);
  std::for_each(vertices_graph.first,
                vertices_graph.second,
                [&assignment]
                (const VertexType& v) {
                  assignment[v] = PortAssignment::UnAssigned;
                });

  // Init open list with a root node
  open_list_t openlist;
  const auto root = *vertices_graph.first;
  openlist.push(root);
  assignment[root] = PortAssignment::PortA;

  while (!openlist.empty()) {
    const auto& node = openlist.front();
    const auto& this_assignment = assignment.at(node);
    auto adjacent_vs = adjacent_vertices(node, spanning_tree);
    std::for_each(adjacent_vs.first,
                  adjacent_vs.second,
                  [&this_assignment, &assignment, &openlist]
                  (const VertexType& adj_v) {
                    auto& adj_assignment = assignment.at(adj_v);

                    switch (adj_assignment) {
                      case PortAssignment::UnAssigned:
                        if (this_assignment == PortAssignment::PortA) {
                          adj_assignment = PortAssignment::PortB;
                        } else {
                          adj_assignment = PortAssignment::PortA;
                        }
                        openlist.push(adj_v);
                        break;

                      case PortAssignment::PortA:
                      case PortAssignment::PortB:
                        if (adj_assignment == this_assignment) {
                          throw std::runtime_error(
                              "Malformed spanning tree!");
                        }
                        break;

                      default:
                        throw std::runtime_error(
                            "Impossibile configuration in the assignment");
                    }
                  });
    openlist.pop();
  }
}

}  // namespace pap_solver


#endif  // __PAP_BCO_PARSER__ALGORITHM__HPP
