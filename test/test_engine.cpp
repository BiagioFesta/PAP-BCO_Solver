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

#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include "../src/matrix_parser.hpp"
#include "../src/Engine.hpp"
using pap_solver::MatrixParser;
using pap_solver::Engine;

typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS> Graph;

static const std::string embedded_matrix =
    "111\n"
    "11\n"
    "1";

template<typename Graph>
static void parse_matrix(Graph* graph) {
  static MatrixParser parser;
  std::stringstream stream;
  stream.str(embedded_matrix);

  auto add_edge_function = [](const typename Graph::vertex_descriptor& v1,
                              const typename Graph::vertex_descriptor& v2,
                              Graph* g) {
    boost::add_edge(v1, v2, *g);
  };
  parser.parse_compressed_matrix(&stream, graph, add_edge_function);
  std::cout << "Number of vertices: " << boost::num_vertices(*graph) << '\n';
  std::cout << "Number of edges: " << boost::num_edges(*graph) << '\n';
}

int main(int argc, char *argv[]) {
  Graph graph;
  Engine<Graph> engine;
  parse_matrix(&graph);
  engine.find_a_solution_and_print(graph, &std::cout);

  return 0;
}
