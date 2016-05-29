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
#include <istream>
#include <boost/graph/adjacency_list.hpp>
#include "options.hpp"

namespace pap_solver {

class PAP_BCO_Solver {
 public:
  ///! A Graph type.
  typedef boost::adjacency_list<boost::vecS,
                                boost::vecS,
                                boost::undirectedS> Graph;
  ///! Default constructor
  PAP_BCO_Solver(int argc, char* argv[]);

  void run();

 private:
  struct AddEdge {
    void operator()(size_t v1, size_t v2, Graph* g) {
      boost::add_edge(v1, v2, *g);
    }
  };
  ProgramOptions m_options;
  Graph m_graph;

  ///! \brief Parse the program option through the command line.
  int parse_cmdline_options(int argc, char* argv[]);

  ///! \brief Parse the m_graph reading the input file.
  void parse_matrix_fromfile();
};

}  // namespace pap_solver

#endif
