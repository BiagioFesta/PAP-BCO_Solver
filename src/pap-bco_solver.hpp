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
#include "GraphUtility.hpp"

namespace pap_solver {

class PAP_BCO_Solver {
 public:
  /// A Graph type.
  typedef boost::adjacency_list<boost::vecS,
                                boost::vecS,
                                boost::undirectedS> Graph;
  /// Default constructor
  PAP_BCO_Solver() noexcept;

  void run(int argc, char* argv[]);

 private:
  ProgramOptions m_options;
  Graph m_graph;

  /// @brief Parse the program option through the command line.
  int parse_cmdline_options(int argc, char* argv[]);

  /// @brief Parse the m_graph reading the input file or stdin.
  void parse_matrix();

  /// @brief Generate a random matrix and prints it
  ///        on the std output or on the file if it has been
  ///        specified in the options.
  template<typename RND>
  void generate_random_matrix(RND* rnd_engine) const;

  /// @brief Prints a briefly help guide on the std output.
  void print_help() const noexcept;

  /// @brief It prints a briefly header on the std output.
  void print_header() const noexcept;
};


template<typename RND>
void PAP_BCO_Solver::generate_random_matrix(RND* rnd_engine) const {
  std::ostream* os = &std::cout;
  std::ofstream file;
  if (m_options.input_filename.size() != 0) {
    file.open(m_options.input_filename);
    os = &file;
  }

  Graph local_graph;
  GraphUtility::generate_random_graph(m_options.generate_num_vertices,
                                      m_options.generate_num_edges,
                                      &local_graph,
                                      rnd_engine);

  GraphUtility::printGraph_asArchive(local_graph, os);

  if (m_options.input_filename.size() != 0) {
    file.close();
  }
}

}  // namespace pap_solver

#endif
