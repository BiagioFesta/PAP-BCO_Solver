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

  /// @brief An map assignment.
  typedef std::map<VertexType, PortAssignment> MapAssignment;

  /// @brief Spanning Tree of Graph type.
  typedef SpanningTree<Graph> SpanningTreeT;

  /// @brief A map<EdgeType, bool> filter.
  typedef typename SpanningTree<Graph>::EdgeFilter EdgeFilter;

  /// @brief A predicate for EdgeFilter
  typedef typename SpanningTree<Graph>::PredicateFilterEdge PredicateFilterEdge;

  /// @brief A view on a graph.
  typedef typename SpanningTree<Graph>::FilteredGraph FilteredGraph;

  /// @brief A struct with the solution and some its details.
  struct Solution {
    MapAssignment m_assignment;
    size_t m_size_solution;
    SpanningTree<Graph> m_spanning_tree;
    std::chrono::milliseconds m_time_for_solution;
  };

  /// @brief Default constructor.
  Algorithm() = default;

  void solve(const Graph& graph, Solution* out_solution);

  /// @brief Sets the seed for the random engine.
  void set_seed(int seed);

  static void print_solution(std::ostream* os,
                             const Solution& solution);

  /// @brief Compute the cutset for the edges which don't belong
  ///        to the spanning tree.
  ///
  /// @param [in] graph                  The main graph
  /// @param [in]current_solution        The solution found
  /// @param edge_of_spanning_tree       The edge which perform on the
  ///                                    cut set. Must to be in tree.
  /// @param cutset_output               The result.
  ///
  /// @note Complexity should be O(|V| + |E|).
  static void fundamental_cutset(const Graph& graph,
                                 const Solution& current_solution,
                                 const EdgeType& edge_of_spanning_tree,
                                 EdgeFilter* cutset_output);

  void debug(const Graph& graph);

 private:
  RndGenerator m_rnd_engine;

  void solve_problem(const Graph& graph, Solution* out_solution);

  static size_t find_all_odd_cotree_edges(const Graph& graph,
                                          const EdgeFilter& edges_in_tree,
                                          EdgeFilter* odd_edges);

  static bool is_odd_cotree_edge(const Graph& graph,
                                 const EdgeFilter& edges_in_spanning_tree,
                                 const EdgeType& e_to_test);

  static void assign_accordance_odd_cycle(const Graph& graph,
                                          EdgeFilter* odd_cotree_edges,
                                          Solution* output_solution);

  /// TODO(biagio): comment
  ///
  /// @param [in] graph               The main graph.
  /// @param [in,out] out_solution    The solution where the output will
  ///                                 be stored.
  ///
  /// @note The spanning tree in the solution must to be already generated.
  static void assign_accordance_spanning_tree(const Graph& graph,
                                              Solution* out_solution);
};

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::set_seed(int seed) {
  m_rnd_engine.seed(seed);
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::solve(
    const Graph& graph,
    Solution* out_solution) {
  assert(out_solution != nullptr);

  auto time_start = std::chrono::steady_clock::now();
  solve_problem(graph, out_solution);
  auto time_stop = std::chrono::steady_clock::now();

  out_solution->m_time_for_solution =
      std::chrono::duration_cast<std::chrono::milliseconds>(
      time_stop - time_start);
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::solve_problem(
    const Graph& graph, Solution* out_solution) {
  assert(out_solution != nullptr);

  auto& assignment = out_solution->m_assignment;

  // Generate the spanning tree
  out_solution->m_spanning_tree.generate_rnd_spanning_tree(graph,
                                                           &m_rnd_engine);

  // Preliminary assignment in according to the spanning tree
  Algorithm::assign_accordance_spanning_tree(graph, out_solution);

  // Find odd cotree edges
  const auto& edges_in_spanning_tree =
      out_solution->m_spanning_tree.get_edges_in_spanning_tree();
  EdgeFilter odd_cotree_edges;
  find_all_odd_cotree_edges(graph, edges_in_spanning_tree, &odd_cotree_edges);

  // Assign V_ab group in according to odd cycles
  assign_accordance_odd_cycle(graph, &odd_cotree_edges, out_solution);
}

template<typename Graph, typename RndGenerator>
bool Algorithm<Graph, RndGenerator>::is_odd_cotree_edge(
    const Graph& graph, const EdgeFilter& edges_in_spanning_tree,
    const EdgeType& e_to_test) {

  auto edges_in_spanning_tree_plus_e = edges_in_spanning_tree;
  edges_in_spanning_tree_plus_e[e_to_test] = true;

  FilteredGraph cycled_graph(
      graph, PredicateFilterEdge(edges_in_spanning_tree_plus_e));

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
void Algorithm<Graph, RndGenerator>::assign_accordance_spanning_tree(
    const Graph& graph, Solution* out_solution) {
  typedef std::queue<VertexType> open_list_t;

  auto& assignment = out_solution->m_assignment;
  const auto& spanning_tree =
      out_solution->m_spanning_tree.get_filtered_graph(graph);

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

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::print_solution(
    std::ostream* os, const Solution& solution) {
  assert(os != nullptr);

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
void Algorithm<Graph, RndGenerator>::fundamental_cutset(
    const Graph& graph, const Solution& current_solution,
    const EdgeType& edge_of_spanning_tree, EdgeFilter* cutset_output) {
  assert(cutset_output != nullptr);

  // Types definitions
  typedef boost::two_bit_color_map<> color_map_t;
  typedef std::queue<VertexType> open_list_t;

  // Colors definitions
  static const auto white_t =  boost::color_traits<
    boost::two_bit_color_type>::white();
  static const auto black_t =  boost::color_traits<
    boost::two_bit_color_type>::black();
  static const auto gray_t =  boost::color_traits<
    boost::two_bit_color_type>::gray();

  const auto& edges_tree = current_solution.
      m_spanning_tree.get_edges_in_spanning_tree();

  // Check if edge is in spanning_tree
  // TODO(biagio): in release puoi levare questo check?
  if (edges_tree.at(edge_of_spanning_tree) == false) {
    throw std::runtime_error("Try to perform cutset of co-tree edge!");
  }

  EdgeFilter& edges_in_cutset = *cutset_output;
  edges_in_cutset.clear();

  // Vertices of edge we're considerating
  auto target = boost::target(edge_of_spanning_tree, graph);
  auto source = boost::source(edge_of_spanning_tree, graph);

  // Some variables inizialization
  const auto num_vertices = boost::num_vertices(graph);
  color_map_t color_map(num_vertices);
  open_list_t openlist;
  const auto spanning_tree_graph = current_solution.
      m_spanning_tree.get_filtered_graph(graph);

  // TODO(biagio): le due colorazioni potrebbero essere parallelizzate
  // Black coloration
  put(color_map, target, black_t);
  openlist.push(target);
  while (!openlist.empty()) {
    const auto& node = openlist.front();
    auto adj_vertices = adjacent_vertices(node, spanning_tree_graph);
    std::for_each(adj_vertices.first,
                  adj_vertices.second,
                  [&source, &color_map, &openlist]
                  (const VertexType& v) {
                    if (v != source) {
                      const auto& v_color = get(color_map, v);
                      if (v_color == white_t) {
                        put(color_map, v, black_t);
                        openlist.push(v);
                      }
                      assert(v_color != gray_t);
                    }
                  });
    openlist.pop();
  }


  // Gray coloration
  put(color_map, source, gray_t);
  openlist.push(source);
  while (!openlist.empty()) {
    const auto& node = openlist.front();
    auto adj_vertices = boost::adjacent_vertices(node, spanning_tree_graph);
    std::for_each(adj_vertices.first,
                  adj_vertices.second,
                  [&target, &color_map, &openlist]
                  (const VertexType& v) {
                    if (v != target) {
                      const auto& v_color = get(color_map, v);
                      if (v_color == white_t) {
                        put(color_map, v, gray_t);
                        openlist.push(v);
                      }
                      assert(v_color != black_t);
                    }
                  });
    openlist.pop();
  }

  // Cut assignment
  const auto edges_spanning = edges(graph);
  std::for_each(edges_spanning.first,
                edges_spanning.second,
                [&edges_in_cutset, &color_map, &graph]
                (const EdgeType& e) {
                  const auto source = boost::source(e, graph);
                  const auto target = boost::target(e, graph);
                  const auto color_s = get(color_map, source);
                  const auto color_t = get(color_map, target);
                  if (color_s != color_t) {
                    edges_in_cutset[e] = true;
                  } else {
                    edges_in_cutset[e]= false;
                  }
                });
}

template<typename Graph, typename RndGenerator>
size_t Algorithm<Graph, RndGenerator>::find_all_odd_cotree_edges(
    const Graph& graph,
    const EdgeFilter& edges_in_tree,
    EdgeFilter* odd_edges) {
  assert(odd_edges != nullptr);

  EdgeFilter& odd_cotree_edges = *odd_edges;
  odd_cotree_edges.clear();
  size_t number_of_odd = 0;
  std::for_each(edges(graph).first,
                edges(graph).second,
                [&graph, &odd_cotree_edges, &edges_in_tree, &number_of_odd]
                (const EdgeType& e) {
                  const auto finder = edges_in_tree.find(e);
                  if (edges_in_tree.at(e) == false) {
                    const auto is_odd =
                        Algorithm::is_odd_cotree_edge(graph,
                                                      edges_in_tree,
                                                      e);
                    if (is_odd == true) {
                      ++number_of_odd;
                    }
                    odd_cotree_edges[e] = is_odd;
                  } else {
                    odd_cotree_edges[e] = false;
                  }
                });
  return number_of_odd;
}


template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::debug(const Graph& graph) {
  set_seed(std::chrono::system_clock::now().time_since_epoch().count());
  Solution solution;

  solution.m_spanning_tree.generate_rnd_spanning_tree(graph, &m_rnd_engine);
  //solution.m_spanning_tree.print(graph, &std::cout);
  const EdgeFilter& edges_tree =
      solution.m_spanning_tree.get_edges_in_spanning_tree();

  EdgeFilter odd_edges;
  auto num_odd = find_all_odd_cotree_edges(graph, edges_tree, &odd_edges);
  std::cout << "Number odd: " << num_odd << '\n';
  
  EdgeFilter cutset;
  const EdgeType* max_current = nullptr;
  int max_odd_cutset = 0;
  EdgeFilter max_cutset;
  for (const auto& e : edges_tree) {
    if (e.second == true) {
      fundamental_cutset(graph, solution, e.first, &cutset);
      int odd_in_cutset = 0;
      std::for_each(cutset.cbegin(),
                    cutset.cend(),
                    [&odd_in_cutset, &odd_edges]
                    (const typename EdgeFilter::value_type& pair) {
                      if (pair.second == true) {
                        if (odd_edges.at(pair.first) == true) {
                          ++odd_in_cutset;
                        } else {
                          --odd_in_cutset;
                        }
                      }
                    });
      //std::cout << e.first << "  :  ";
      //std::cout << odd_in_cutset << '\n';
      if (max_current == nullptr) {
        max_current = &(e.first);
        max_odd_cutset = odd_in_cutset;
        max_cutset = cutset;
      } else {
        if (odd_in_cutset > max_odd_cutset) {
          max_current = &(e.first);
          max_odd_cutset = odd_in_cutset;
          max_cutset = cutset;
        }
      }
    }
  }
  std::cout << "The max is: " << *max_current << '\n';
  assign_accordance_spanning_tree(graph, &solution);
  assign_accordance_odd_cycle(graph, &odd_edges, &solution);
  std::cout << "Size of solution: " << solution.m_size_solution << '\n';
  
  solution.m_spanning_tree.perform_transformation(
      graph,
      std::find_if(max_cutset.cbegin(),
                   max_cutset.cend(),
                   [&odd_edges]
                   (const typename decltype(max_cutset)::value_type& p) {
                     return p.second && odd_edges.at(p.first) == true;
                   })->first,
      *max_current);
  const auto& edges_tree2 =
      solution.m_spanning_tree.get_edges_in_spanning_tree();
  //solution.m_spanning_tree.print(graph, &std::cout);
  num_odd = find_all_odd_cotree_edges(graph, edges_tree2, &odd_edges);
  std::cout << "Number odd: " << num_odd << '\n';
  assign_accordance_spanning_tree(graph, &solution);
  assign_accordance_odd_cycle(graph, &odd_edges, &solution);
  std::cout << "Size solution: " << solution.m_size_solution << '\n';
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::assign_accordance_odd_cycle(
    const Graph& graph,
    EdgeFilter* oddCotreeEdges,
    Solution* output_solution) {
  assert(output_solution != nullptr);

  auto& assignment = output_solution->m_assignment;
  auto odd_cotree_edges = *oddCotreeEdges;

  // Create a g' considering only odd cotree edges
  FilteredGraph g_prime(graph, PredicateFilterEdge(odd_cotree_edges));

  // Solve odd edges's problem
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
        assignment[target] = PortAssignment::PortAB;
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
        assignment[max_vertex_degree] = PortAssignment::PortAB;
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
  output_solution->m_size_solution = number_of_AB;
}


}  // namespace pap_solver


#endif  // __PAP_BCO_PARSER__ALGORITHM__HPP
