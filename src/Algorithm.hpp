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
#include <limits>
#include <queue>
#include <ostream>
#include <utility>
#include <chrono>
#include <stack>
#include <vector>
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

  /// @brief An elementary transformation on a spanning tree.
  typedef std::pair<EdgeType, EdgeType> EdgeTransformation;

  /// @brief A struct with the solution and some its details.
  struct Solution {
    MapAssignment m_assignment;
    size_t m_size_solution;
    SpanningTree<Graph> m_spanning_tree;
    size_t m_number_odd_edges;
    std::chrono::milliseconds m_time_for_solution;
  };

  /// @brief Default constructor.
  Algorithm() = default;

  void find_rnd_solution_fast(const Graph& graph, Solution* out_solution);

  void best_local_solution(const Graph& graph, Solution* out_solution);

  void best_local_solution_h(const Graph& graph, Solution* out_solution);

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
                                 const SpanningTreeT& spanning_tree,
                                 const EdgeType& edge_of_spanning_tree,
                                 EdgeFilter* cutset_output);

 private:
  struct InstanceSolution {
    MapAssignment m_assignment;
    size_t m_size_solution;
    EdgeFilter m_odd_edges;
    size_t m_num_odd_edges;

    /// @brief This is an intermediate solution based on
    /// assignment of vertices respect with the spanning tree.
    /// This assignment is usefull because is used by the
    /// method to discover odd edges.
    MapAssignment m_mapped_sp_based;
  };

  RndGenerator m_rnd_engine;

  static void solve_problem_for_a_tree(const Graph& graph,
                                       const SpanningTreeT& spanning_tree,
                                       InstanceSolution* output_solution);

  static size_t find_all_odd_cotree_edges(const Graph& graph,
                                          const MapAssignment& mapped_sp_based,
                                          EdgeFilter* odd_edges);

  static bool is_odd_cotree_edge(const Graph& graph,
                                 const MapAssignment& mapped_sp_based,
                                 const EdgeType& e_to_test);

  /// @brief Assign some vertices to the portAB to solve_problem the problem
  ///        of solve_problem the problem of odd cycles.
  ///
  /// @param [in] graph                   The principal graph.
  /// @param [in,out] odd_cotree_edges    A filter of all odd cotree edges.
  ///                                     This filter will be destroy by the
  ///                                     function.
  /// @param [out] assignment             A valid assignment.
  /// @param [out] number_of_AB           The number of assignment AB
  ///                                     in the result of the function.
  static void assign_accordance_odd_cycle(const Graph& graph,
                                          EdgeFilter* odd_cotree_edges,
                                          MapAssignment* assignment,
                                          size_t* number_of_AB);

  /// @brief Assign all vertices to the port A or the port B in according
  ///        to the spanning tree passed as argument.
  ///
  /// @param [in] graph             The principal graph.
  /// @param spanning_tree_graph    The spanning tree of graph.
  /// @param assignment             The assignment in output.
  static void assign_accordance_spanning_tree(
      const Graph& graph,
      const FilteredGraph& spanning_tree_graph,
      MapAssignment* assignment);

  /// @brief The function find the best transformation for a random
  ///        spanning tree.
  /// @param [in] graph                The main graph of the problem
  /// @param st                        A random spanning tree.
  ///                                  Note that you don't have to generate
  ///                                  the solution for that tree.
  /// @param [out] out                 The best local solution for that tree.
  /// @param [out] out_transformation  The best transformation you can perform
  ///                                  on the tree in order to achive the
  ///                                  best local solution.
  /// @return 'true' if a transformation has been found, and it will be
  ///    stored in the param 'out_transformation', otherwise 'false'
  ///    that will mean the out param will be useless.
  static bool find_best_local_solution_inTree(
      const Graph& graph,
      const SpanningTreeT& st,
      InstanceSolution* out,
      EdgeTransformation* out_transformation);

  static bool find_best_local_solution_inTree_heristic(
      const Graph& graph,
      const SpanningTreeT& st,
      InstanceSolution* out,
      EdgeTransformation* out_transformation);

  /// @brief Cycle detection algorithm for UNDIRECTED GRAPH!
  ///
  /// @param [in] graph          An undirected graph. VertexType must to be
  ///                            an integer type.
  /// @param [in] num_vertices   The number of vertices in the graph.
  /// @param [out] edge_in_cycle An EdgeDescriptor which will be an
  ///                            edge in the cycle.
  ///
  /// @return true whether there is a cycle in the graph.
  /// @note The output parameter will be modified even if the result will
  /// be false. In that case the value will be no sense.
  template<typename GenericGraph, typename EdgeDescriptor>
  static bool a_cycle_detection(const GenericGraph& graph,
                                const size_t num_vertices,
                                EdgeDescriptor* edge_in_cycle);
};

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::set_seed(int seed) {
  m_rnd_engine.seed(seed);
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::find_rnd_solution_fast(
    const Graph& graph, Solution* out_solution) {
  assert(out_solution != nullptr);

  // Local variables
  SpanningTreeT rnd_spanning_tree;
  InstanceSolution local_solution;

  auto time_start = std::chrono::steady_clock::now();

  // Generate a rnd spanning tree
  rnd_spanning_tree.generate_rnd_spanning_tree(graph, &m_rnd_engine);

  // Solve the problem
  solve_problem_for_a_tree(graph,
                           rnd_spanning_tree,
                           &local_solution);

  auto time_stop = std::chrono::steady_clock::now();

  // Assign local solution on the output
  // TODO(biagio): verificare che la move sia implementata e che non copy
  out_solution->m_spanning_tree = std::move(rnd_spanning_tree);
  out_solution->m_size_solution = std::move(local_solution.m_size_solution);
  out_solution->m_assignment = std::move(local_solution.m_assignment);
  out_solution->m_number_odd_edges = std::move(local_solution.m_num_odd_edges);

  // Assign the time elapsed
  out_solution->m_time_for_solution =
      std::chrono::duration_cast<std::chrono::milliseconds>(time_stop -
                                                            time_start);
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::best_local_solution(
    const Graph& graph, Solution* out_solution) {
  assert(out_solution != nullptr);

  // Function setting
  static constexpr size_t NUMBER_OF_TREE_TO_GENERATE = 10;

  // Local variables
  SpanningTreeT rnd_spanning_tree;
  InstanceSolution local_solution;
  EdgeTransformation local_transformation;
  out_solution->m_size_solution = std::numeric_limits<size_t>::max();

  auto time_start = std::chrono::steady_clock::now();

  for (size_t i = 0; i < NUMBER_OF_TREE_TO_GENERATE; ++i) {
    // Generate a rnd spanning tree
    rnd_spanning_tree.generate_rnd_spanning_tree(graph, &m_rnd_engine);

    // Rafine the tree at the best you can
    while (find_best_local_solution_inTree_heristic(graph, rnd_spanning_tree,
                                                    &local_solution,
                                                    &local_transformation) == true) {
      rnd_spanning_tree.perform_transformation(graph,
                                               local_transformation.first,
                                               local_transformation.second);
      if (local_solution.m_size_solution < out_solution->m_size_solution) {
        out_solution->m_size_solution = local_solution.m_size_solution;
        out_solution->m_assignment = std::move(local_solution.m_assignment);
        out_solution->m_number_odd_edges = local_solution.m_num_odd_edges;
      }
    }

    if (local_solution.m_size_solution < out_solution->m_size_solution) {
      out_solution->m_size_solution = local_solution.m_size_solution;
      out_solution->m_assignment = std::move(local_solution.m_assignment);
      out_solution->m_number_odd_edges = local_solution.m_num_odd_edges;
    }
  }

  auto time_stop = std::chrono::steady_clock::now();
  // Assign the time elapsed
  out_solution->m_time_for_solution =
      std::chrono::duration_cast<std::chrono::milliseconds>(time_stop -
                                                            time_start);
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::solve_problem_for_a_tree(
    const Graph& graph,
    const SpanningTreeT& spanning_tree,
    InstanceSolution* output_solution) {
  assert(output_solution != nullptr);

  const auto& spanning_tree_graph =
      FilteredGraph(graph,
                    PredicateFilterEdge(spanning_tree.
                                        get_edges_in_spanning_tree()));
  const auto& edges_in_spanning_tree =
      spanning_tree.get_edges_in_spanning_tree();

  // Preliminary assignment in according to the spanning tree
  assign_accordance_spanning_tree(graph,
                                  spanning_tree_graph,
                                  &output_solution->m_assignment);

  // Fill mapped_spanning_tree-based
  output_solution->m_mapped_sp_based = output_solution->m_assignment;

  // Find odd cotree edges
  output_solution->m_num_odd_edges =
      find_all_odd_cotree_edges(graph,
                                output_solution->m_assignment,
                                &output_solution->m_odd_edges);

  // Assign V_ab group in according to odd cycles
  EdgeFilter copy_odd_edges = output_solution->m_odd_edges;
  assign_accordance_odd_cycle(graph,
                              &copy_odd_edges,
                              &output_solution->m_assignment,
                              &output_solution->m_size_solution);
}

template<typename Graph, typename RndGenerator>
bool Algorithm<Graph, RndGenerator>::is_odd_cotree_edge(
    const Graph& graph,
    const MapAssignment& mapped_sp_based,
    const EdgeType& e_to_test) {
  const auto& target = boost::target(e_to_test, graph);
  const auto& source = boost::source(e_to_test, graph);
  assert(mapped_sp_based.at(target) != PortAssignment::UnAssigned);
  assert(mapped_sp_based.at(source) != PortAssignment::UnAssigned);

  return mapped_sp_based.at(target) == mapped_sp_based.at(source) ?
      true : false;
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::assign_accordance_spanning_tree(
    const Graph& graph,
    const FilteredGraph& spanning_tree_graph,
    MapAssignment* p_assignment) {
  assert(p_assignment != nullptr);

  typedef std::queue<VertexType> open_list_t;

  auto& assignment = *p_assignment;

  // Initialization all vertices are unassigned O(V)
  const auto vertices_graph = vertices(spanning_tree_graph);
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
    auto adjacent_vs = adjacent_vertices(node, spanning_tree_graph);
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
    const Graph& graph, const SpanningTreeT& spanning_tree,
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

  const auto& edges_tree = spanning_tree.get_edges_in_spanning_tree();

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
  const auto spanning_tree_graph = spanning_tree.get_filtered_graph(graph);

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
  const auto edges_graph = edges(graph);
  std::for_each(edges_graph.first,
                edges_graph.second,
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
    const MapAssignment& mapped_sp_based,
    EdgeFilter* odd_edges) {
  assert(odd_edges != nullptr);

  EdgeFilter& odd_cotree_edges = *odd_edges;
  odd_cotree_edges.clear();
  size_t number_of_odd = 0;

  bool local_odd;
  const auto edges_graph = edges(graph);
  std::for_each(edges_graph.first,
                edges_graph.second,
                [&]
                (const EdgeType& e) {
                 local_odd =
                     Algorithm::is_odd_cotree_edge(graph,
                                                   mapped_sp_based,
                                                   e);
                 if (local_odd == true) {
                   ++number_of_odd;
                 }

                 odd_cotree_edges[e] = local_odd;
                });
  return number_of_odd;
}

template<typename Graph, typename RndGenerator>
void Algorithm<Graph, RndGenerator>::assign_accordance_odd_cycle(
    const Graph& graph,
    EdgeFilter* p_odd_cotree_edges,
    MapAssignment* p_assignment,
    size_t* p_number_of_AB) {
  assert(p_odd_cotree_edges != nullptr);
  assert(p_assignment != nullptr);
  assert(p_number_of_AB != nullptr);

  auto& assignment = *p_assignment;
  auto& odd_cotree_edges = *p_odd_cotree_edges;
  auto& number_of_AB = *p_number_of_AB;

  // Create a g' considering only odd cotree edges
  FilteredGraph g_prime(graph, PredicateFilterEdge(odd_cotree_edges));

  // Solve odd edges's problem
  bool found_one_degree;
  bool g_prime_empty = false;
  number_of_AB = 0;
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
}

template<typename Graph, typename RndGenerator>
bool Algorithm<Graph, RndGenerator>::find_best_local_solution_inTree(
    const Graph& graph, const SpanningTreeT& st,
    InstanceSolution* out, EdgeTransformation* out_transformation) {
  assert(out != nullptr);
  assert(out_transformation != nullptr);

  // Local variables
  EdgeFilter cutset;
  SpanningTreeT sp_copy = st;
  InstanceSolution local_solution;
  bool rts =  false;

  // Find the initial solution
  solve_problem_for_a_tree(graph, st, out);

  // View on spanning tree
  const auto&& sp_graph =
      FilteredGraph(graph,
                    PredicateFilterEdge(st.get_edges_in_spanning_tree()));

  // For each edge in the spanning tree
  const auto& all_edges = edges(sp_graph);
  std::for_each(all_edges.first,
                all_edges.second,
                [&]
                (const EdgeType& e_sp) {
                  fundamental_cutset(graph, st, e_sp, &cutset);

                  // For each edge in the cuset
                  std::for_each(
                      cutset.cbegin(),
                      cutset.cend(),
                      [&]
                      (const typename EdgeFilter::value_type& pe_cuset) {
                        // This edge of the cutset
                        const auto& e_ct = pe_cuset.first;

                        // Check if the edge belong to the cutset and
                        // it is odd
                        if (pe_cuset.second == true &&
                            out->m_odd_edges[e_ct] == true &&
                            e_ct != e_sp) {
                          // Perform the transformation
                          sp_copy.perform_transformation(graph, e_ct, e_sp);

                          // Solve the problem for the new tree
                          solve_problem_for_a_tree(graph,
                                                   sp_copy,
                                                   &local_solution);

                          // Check whether the new solution is better
                          if (local_solution.m_size_solution <
                              out->m_size_solution) {
                            // Set the new solution
                            out->m_size_solution =
                                local_solution.m_size_solution;
                            out->m_assignment = std::move(
                                local_solution.m_assignment);
                            out->m_odd_edges = std::move(
                                local_solution.m_odd_edges);
                            out->m_num_odd_edges =
                                local_solution.m_num_odd_edges;
                            out->m_mapped_sp_based = std::move(
                                local_solution.m_mapped_sp_based);

                            // Set the transformation as better so far
                            *out_transformation = std::make_pair(e_ct, e_sp);

                            // Change the flag result
                            rts = true;
                          }

                          // Reverse the transformation
                          sp_copy.perform_transformation(graph, e_sp, e_ct);
                        }
                      });
                });
  return rts;
}

template<typename Graph, typename RndGenerator>
template<typename GenericGraph, typename EdgeDescriptor>
bool Algorithm<Graph, RndGenerator>::a_cycle_detection(
    const GenericGraph& graph,
    const size_t num_vertices,
    EdgeDescriptor* edge_in_cycle) {
  assert(edge_in_cycle != nullptr);
  if (num_vertices == 0) {
    throw std::invalid_argument(
        "Cycle detection algorithm with empty graph has been called.");
  }

  // Type declaration
  typedef typename GenericGraph::vertex_descriptor Vertex;
  static_assert(std::is_integral<Vertex>::value,
                "The vertex type (descriptor) must to be a integer type!");
  typedef std::stack<std::pair<Vertex, Vertex>> OpenList;
  // TODO(biagio): essendo integere i vertici puoi anche usare un array
  // indicizzato dai vertici stessi! Lo spazio sarebbe lo stesso
  typedef std::unordered_map<Vertex, void*> ClosedList;

  // Algorithm idea: Depth First Search for at most the number of vertices.
  // Since at most n − 1 edges can be tree edges, a cycle has to be
  // discovered before UNDIRECTED GRAPH!

  OpenList open_list;
  ClosedList closed_list;
  closed_list.reserve(num_vertices);

  const auto& vertices = boost::vertices(graph);
  open_list.push(std::make_pair(*vertices.first, *vertices.second));

  bool cycle_detected = false;
  while (open_list.empty() == false && cycle_detected == false) {
    const auto current_node = open_list.top().first;
    const auto current_parent = open_list.top().second;
    closed_list[current_node] = nullptr;
    open_list.pop();

    const auto adj_vertices = boost::adjacent_vertices(current_node, graph);
    std::for_each(adj_vertices.first,
                  adj_vertices.second,
                  [&] (const Vertex& adj_v) {
                    if (closed_list.find(adj_v) == closed_list.cend()) {
                      open_list.push(std::make_pair(adj_v, current_node));
                    } else if (adj_v != current_parent) {
                      auto edge_link =
                          boost::edge(adj_v, current_node, graph);
                      assert(edge_link.second == true);
                      *edge_in_cycle = edge_link.first;
                      cycle_detected = true;
                    }
                  });

    // Since the graph could be disconnected
    if (cycle_detected == false &&
        open_list.empty() == true &&
        closed_list.size() < num_vertices) {
      // No cycle but there still are vertices to explore
      const auto finder = std::find_if(
          vertices.first, vertices.second,
          [&] (const Vertex& v) {
            if (closed_list.find(v) == closed_list.cend()) {
              return true;
            }
            return false;
          });
      assert(finder != vertices.second);
      open_list.push(std::make_pair(*finder, *vertices.second));
    }
  }

  return cycle_detected;
}

template<typename Graph, typename RndGenerator>
bool Algorithm<Graph, RndGenerator>::find_best_local_solution_inTree_heristic(
    const Graph& graph,
    const SpanningTreeT& st,
    InstanceSolution* out,
    EdgeTransformation* out_transformation) {
  assert(out != nullptr);
  assert(out_transformation != nullptr);

  // Local variables
  EdgeFilter cutset;
  EdgeType local_edge;
  SpanningTreeT sp_copy = st;
  InstanceSolution local_solution;
  const auto num_vertices = boost::num_vertices(graph);
  const auto num_edges = boost::num_edges(graph);

  // Find the initial solution
  solve_problem_for_a_tree(graph, st, out);

  // View on spanning tree
  const auto& edges_in_tree = st.get_edges_in_spanning_tree();
  const auto&& sp_graph =
      FilteredGraph(graph,
                    PredicateFilterEdge(edges_in_tree));

  // View on odd co-tree graph
  const auto&& odd_graph =
      FilteredGraph(graph,
                    PredicateFilterEdge(out->m_odd_edges));

  // Find a cycle in the odd graph
  if (a_cycle_detection(odd_graph, num_vertices, &local_edge)) {
    // Find all edges we could swap
    std::vector<EdgeType> cand_edges_to_swap;

    // TODO(biagio): troppo pessimistico!
    cand_edges_to_swap.reserve(num_edges);

    // TODO(biagio): this algorithm is a brute force. No good performance
    auto edges_in_stree = boost::edges(sp_graph);
    std::for_each(
        edges_in_stree.first, edges_in_stree.second,
        [&] (const EdgeType& e) {
          fundamental_cutset(graph, st, e, &cutset);

          const auto finder = std::find_if(
              cutset.cbegin(), cutset.cend(),
              [&] (const typename decltype(cutset)::value_type& p) {
                if (p.second == true && p.first == local_edge) {
                  return true;
                }
                return false;
              });

          if (finder != cutset.cend()) {
            cand_edges_to_swap.push_back(e);
          }
        });

    for (const auto& e : cand_edges_to_swap) {
      // Apply transformation in order to try to remove the cycle.
      sp_copy.perform_transformation(graph, local_edge, e);

      // Solve the problem for the new tree
      solve_problem_for_a_tree(graph,
                               sp_copy,
                               &local_solution);

      // Check whether the new solution is better
      if (local_solution.m_size_solution < out->m_size_solution) {
        // Set the new solution
        out->m_size_solution =
            local_solution.m_size_solution;
        out->m_assignment = std::move(
            local_solution.m_assignment);
        out->m_odd_edges = std::move(
            local_solution.m_odd_edges);
        out->m_num_odd_edges =
            local_solution.m_num_odd_edges;
        out->m_mapped_sp_based = std::move(
            local_solution.m_mapped_sp_based);

        // Set the transformation as better so far
        *out_transformation = std::make_pair(local_edge, e);
      }

      // Reverse the transformation
      sp_copy.perform_transformation(graph, e, local_edge);
    }
    return true;
  }
  return false;
}


}  // namespace pap_solver


#endif  // __PAP_BCO_PARSER__ALGORITHM__HPP
