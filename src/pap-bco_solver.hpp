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

#ifndef __PAP_BCO_SOLVER__HPP
#define __PAP_BCO_SOLVER__HPP

#include <string>
#include <random>
#include <boost/graph/adjacency_list.hpp>
#include "options.hpp"
#include "spanning_tree.hpp"
#include "matrix_parser.hpp"

namespace pap_solver {

class PAP_BCO_Solver {
 public:
  /// Properties for a vertex.
  struct VertexProperties {
    enum Port {
      UnAssigned,
      PortA,
      PortB,
      PortAB
    } m_port = Port::UnAssigned;
    bool m_inCoverSet = false;
  };

  /// Properties for a edge.
  struct EdgeProperties {
    bool m_intree = false;
    bool m_odd = false;
  };

  /// A Graph type.
  typedef boost::adjacency_list<boost::vecS,
                                boost::vecS,
                                boost::undirectedS,
                                VertexProperties,
                                EdgeProperties> Graph;
  /// Default constructor
  PAP_BCO_Solver() noexcept;

  void run(int argc, char* argv[]);

 private:
  ProgramOptions m_options;
  Graph m_graph;
  MatrixParser m_mat_parser;

  /// @brief Parse the program option through the command line.
  int parse_cmdline_options(int argc, char* argv[]);

  /// @brief Parse the m_graph reading the input file.
  void parse_matrix_fromfile();

  /// @brief Generate a random matrix and prints it
  ///        on the std output or on the file if it has been
  ///        specified in the options.
  template<typename RND>
  void generate_random_matrix(RND* rnd_engine) const;

  /// @brief Prints a briefly help guide on the std output.
  void print_help() const noexcept;

  /// @brief It prints a briefly header on the std output.
  void print_header() const noexcept;

  void print_all_vertices_and_ports(std::ostream* os) const noexcept;

  void algorithm_assign_port_byTree(const SpanningTree<Graph>& st);

  /// @brief Takes a valid spanning tree and assigns property
  ///        of all adges of the graph. (Such as in_tree, or odd).
  /// @param [in] st    A Valid spanning tree for the graph.
  ///
  void assign_edges_property_byTree(const SpanningTree<Graph>& st);

  ///
  /// @param [in] e    A co-tree edge you want to check.
  /// @param [in] st   A spanning tree of the graph.
  ///
  /// @return whether the co-tree edge ('e') is a odd co-tree edge
  ///         or not.
  /// @note The edge 'e' must to be a co-tree edge.
  bool is_odd_cotree_edge(const Graph::edge_descriptor& e,
                          const SpanningTree<Graph>& st) const;
};

}  // namespace pap_solver

#endif
