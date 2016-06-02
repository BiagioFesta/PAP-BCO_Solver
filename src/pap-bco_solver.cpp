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

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <ctime>
#include <queue>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/bipartite.hpp>
#include "pap-bco_solver.hpp"

int main(int argc, char *argv[]) {
  try {
    pap_solver::PAP_BCO_Solver solver;
    solver.run(argc, argv);
  } catch(const std::exception& err) {
    std::cerr << err.what() << '\n';
    return -1;
  }
  return 0;
}




namespace pap_solver {

PAP_BCO_Solver::PAP_BCO_Solver() noexcept {
}

int PAP_BCO_Solver::parse_cmdline_options(int argc, char* argv[]) {
  static struct option long_options[] = {
    {"compressed", no_argument, 0, 'c'},
    {"help", no_argument, 0, 'h'},
    {"generate", required_argument, 0, 'g'},
    {0, 0, 0, 0}
  };
  int option_index;
  while (true) {
    option_index = 0;
    auto c = getopt_long(argc, argv, "g:hc", long_options, &option_index);
    if (c == -1) break;
    switch (c) {
      case 'c':
        m_options.compressed_matrix = true;
        break;
      case 'h':
        m_options.display_help = true;
        break;
      case 'g':
        m_options.generate_random_matrix = true;
        try {
          m_options.size_generation_matrix = std::stoi(optarg);
        } catch (const std::exception& err) {
          std::cerr <<
              "--generate NUMBER\nNUMBER must to be a integer number!\n";
          return -1;
        }
        break;
      case '?':
        return -1;
      default:
        std::cerr << "Invalid parsing options\n";
        return -1;
    }
  }
  if (optind < argc) {
    m_options.input_filename = argv[optind];
  }
  return 0;
}

void PAP_BCO_Solver::run(int argc, char* argv[]) {
#ifndef _DEBUG
  print_header();
#endif

  if (parse_cmdline_options(argc, argv)) {
    throw std::invalid_argument("Invalid process"
                                " parsing command line arguments");
  }

  if (m_options.display_help == true) {
    print_help();
    return;
  }

  std::default_random_engine rnd_eng(std::time(nullptr));
  // std::default_random_engine rnd_eng(0);

  if (m_options.generate_random_matrix == true) {
    generate_random_matrix(&rnd_eng);
    return;
  }

  if (m_options.input_filename.size() > 0)
    parse_matrix_fromfile();
  else
    parse_matrix_from_stdin();

  SpanningTree<Graph> spanning_tree;
  spanning_tree.makeRandom_fromGraph(m_graph, &rnd_eng);
  spanning_tree.print(&std::cout);
  assign_edges_property_byTree(spanning_tree);
  auto num_AB_vertices = algorithm_assign_port_byTree(spanning_tree);
  print_all_vertices_and_ports(&std::cout);
  std::cout << "Number of portAB: " << num_AB_vertices << '\n';
}

void PAP_BCO_Solver::parse_matrix_fromfile() {
  std::ifstream file;
  file.open(m_options.input_filename);
  if (file.fail() == true) {
    throw std::invalid_argument("The file '" +
                                m_options.input_filename +
                                "' cannot be read.");
  }
  if (m_options.compressed_matrix == true) {
    m_mat_parser.parse_compressed_matrix(&file,
                                         &m_graph,
                                         [](size_t v1, size_t v2, Graph* g) {
                                           boost::add_edge(v1, v2, *g);
                                         });
  } else {
    try {
      m_mat_parser.parse_full_matrix(&file,
                                     &m_graph,
                                     [](size_t v1, size_t v2, Graph* g) {
                                       boost::add_edge(v1, v2, *g);
                                     });
    } catch(const std::exception& err) {
      // TODO(biagio): make the suggestion more usefull
      std::string messages = err.what();
      messages += "\n\tMaybe the option '-c' could be usefull."
          " Use '-h' for help";
      throw std::runtime_error(messages);
    }
  }
  file.close();
}

void PAP_BCO_Solver::parse_matrix_from_stdin() {
  if (m_options.compressed_matrix == true) {
    m_mat_parser.parse_compressed_matrix(&std::cin,
                                         &m_graph,
                                         [](size_t v1, size_t v2, Graph* g) {
                                           boost::add_edge(v1, v2, *g);
                                         });
  } else {
    m_mat_parser.parse_full_matrix(&std::cin,
                                   &m_graph,
                                   [](size_t v1, size_t v2, Graph* g) {
                                     boost::add_edge(v1, v2, *g);
                                   });
  }
}

void PAP_BCO_Solver::print_help() const noexcept {
  std::cout <<
      R"##(
Use: pap-bco_solver [OPTION]... [FILENAME]
   Command Options:
-c, --compresed           : Specify whether the input matrix is a compressed format or not.
-g SIZE, --generate SIZE  : Generate a valid random matrix with dimention SIZE.
-h, --help                : Display this guide.

)##";
}

void PAP_BCO_Solver::print_header() const noexcept {
  std::cout << R"##(
           PAP-BCO Solver
   Copyright (C) 2016  Biagio Festa
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

)##";
}

size_t PAP_BCO_Solver::algorithm_assign_port_byTree(
    const SpanningTree<Graph>& st) {
  const auto& map_tree = st.getMap();

  // First step. Assign each vertex to A or B in according to the spanning
  // tree.
  for (const auto& node : map_tree) {
    auto& port_parent = m_graph[node.second].m_port;
    auto& port_this = m_graph[node.first].m_port;
    if (port_parent == VertexProperties::Port::UnAssigned) {
      port_this = VertexProperties::Port::PortA;
    } else if (port_parent == VertexProperties::Port::PortA) {
      port_this = VertexProperties::Port::PortB;
    } else {
      port_this = VertexProperties::Port::PortA;
    }
  }

  // g' is the ''odd'' co-tree graph. All edges which don't belong to the tree
  // and are odd co-tree edge.
  struct PredEdges {
    PredEdges() = default;
    bool operator()(const Graph::edge_descriptor& e) const {
      // In g' ci sono i nodi che sono co-tree && dispari.
      return (!((*mp_graph)[e].m_intree) && ((*mp_graph)[e].m_odd));
    }
    Graph* mp_graph;
  } pred_edges{&m_graph};
  boost::filtered_graph<Graph, PredEdges> g_prime(m_graph, pred_edges);

  bool found_one_degree;
  bool g_prime_empty = false;
  size_t number_of_AB = 0;

  // Second step of the algorithm.
  while (!g_prime_empty) {
    found_one_degree = false;
    // Looking for an vertex with degree = 1
    for (auto i = boost::vertices(g_prime).first;
         i != boost::vertices(g_prime).second && found_one_degree == false;
         ++i) {
      auto degree = boost::out_degree(*i, g_prime);
      auto i_edges = boost::out_edges(*i, g_prime);
      if (degree == 1) {
        // TODO(biagio): è stupido sai a priori che è uno!
        //              non c'è bisogno for_each
        std::for_each(i_edges.first,
                      i_edges.second,
                      [&g_prime, &number_of_AB]
                      (const Graph::edge_descriptor& e) {
                        auto target = boost::target(e, g_prime);
                        g_prime[target].m_port = VertexProperties::Port::PortAB;
                        ++number_of_AB;
                        g_prime[e].m_intree = true;
                      });
        found_one_degree = true;
      }
    }
    if (found_one_degree == false) {
      // Non ho trovato vertici con degree uguali a 1
      // Quindi cerco quello massimo e lo poto!
      int num_degree_max = 0;
      Graph::vertex_descriptor max_vertex_degree;
      for (auto i = boost::vertices(g_prime).first;
           i != boost::vertices(g_prime).second;
           ++i) {
        auto degree = boost::out_degree(*i, g_prime);
        if (degree > num_degree_max) {
          num_degree_max = degree;
          max_vertex_degree = *i;
        }
      }
      if (num_degree_max > 0) {
        // Ho trovato il massimo, lo aggiunto al coverset e lo  poto!
        g_prime[max_vertex_degree].m_port = VertexProperties::Port::PortAB;
        ++number_of_AB;
        std::for_each(boost::out_edges(max_vertex_degree, g_prime).first,
                      boost::out_edges(max_vertex_degree, g_prime).second,
                      [&g_prime](const Graph::edge_descriptor& e) {
                        g_prime[e].m_intree = true;
                      });
      } else {
        // Non ho trovato massimo => non ci sono edges
        g_prime_empty = true;
      }
    }
  }
  return number_of_AB;
}

void PAP_BCO_Solver::print_all_vertices_and_ports(std::ostream* os)
    const noexcept {
  if (os == nullptr) return;
  auto its = boost::vertices(m_graph);
  for (auto i = its.first; i != its.second; ++i) {
    *os << "Vertex (" << *i << "): ---> ";
    switch (m_graph[*i].m_port) {
      case VertexProperties::Port::PortA:
        *os << "Port A";
        break;
      case VertexProperties::Port::PortB:
        *os << "Port B";
        break;
      case VertexProperties::Port::PortAB:
        *os << "Port AB";
        break;
      case VertexProperties::Port::UnAssigned:
        *os << "Unassigned!";
    }
    *os << '\n';
  }
}

void PAP_BCO_Solver::assign_edges_property_byTree(
    const SpanningTree<Graph>& st) {
  // Set all edges to not belog the spanning tree
  // TODO(biagio): this could be useless because default constructor
  std::for_each(boost::edges(m_graph).first,
                boost::edges(m_graph).second,
                [this](Graph::edge_descriptor e) {
                  m_graph[e].m_intree = false;
                  m_graph[e].m_odd = false;
                });

  // Set in spanning tree property
  for (const auto& node : st.getMap()) {
    std::for_each(boost::out_edges(node.first, m_graph).first,
                  boost::out_edges(node.first, m_graph).second,
                  [this, &node](Graph::edge_descriptor e) {
                    if (boost::target(e, m_graph) == node.second)
                      m_graph[e].m_intree = true;
                  });
  }

  // Set the odd co-tree edge property
  auto its = boost::edges(m_graph);
  for (auto i = its.first; i != its.second; ++i) {
    if (m_graph[*i].m_intree == false) {
      if (is_odd_cotree_edge(*i, st) == true) {
        m_graph[*i].m_odd = true;
      }
    }
  }
}

bool PAP_BCO_Solver::is_odd_cotree_edge(const Graph::edge_descriptor& e,
                                        const SpanningTree<Graph>& st) const {
  typedef boost::two_bit_color_map<> color_map_t;
  typedef std::queue<Graph::vertex_descriptor> open_list_t;
  static const auto white_t =  boost::color_traits<
    boost::two_bit_color_type>::white();
  static const auto black_t =  boost::color_traits<
    boost::two_bit_color_type>::black();
  static const auto gray_t =  boost::color_traits<
    boost::two_bit_color_type>::gray();

  if (m_graph[e].m_intree == true) {
    std::runtime_error("Cannot verify if an edge is odd when it "
                       "belogs the tree");
  }
  const auto num_vertices = boost::num_vertices(m_graph);
  color_map_t color_map(num_vertices);

  open_list_t openlist;
  const auto& root = st.getMap().cbegin()->first;
  openlist.push(root);
  boost::put(color_map, root, black_t);

  while (openlist.empty() == false) {
    const auto& node = openlist.front();
    auto this_node_color = boost::get(color_map, node);
    for (auto i = boost::out_edges(node, m_graph).first;
         i != boost::out_edges(node, m_graph).second;
         ++i) {
      if (m_graph[*i].m_intree == true || *i == e) {
        auto target = boost::target(*i, m_graph);
        auto color_adjacent = boost::get(color_map, target);
        if (color_adjacent == white_t) {
          if (this_node_color == black_t)
            boost::put(color_map, target, gray_t);
          else
            boost::put(color_map, target, black_t);
          openlist.push(target);
        } else if (color_adjacent == this_node_color) {
          return true;
        }
      }
    }
    openlist.pop();
  }
  return false;
}

}  // namespace pap_solver
