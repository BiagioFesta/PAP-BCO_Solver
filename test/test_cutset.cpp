// Copyright 2016 <Biagio Festa>

#include <string>
#include <memory>
#include <cstring>
#include <istream>
#include <fstream>
#include "../src/matrix_parser.hpp"
#include "../src/Engine.hpp"
using pap_solver::MatrixParser;
using pap_solver::Engine;
using pap_solver::SpanningTree;
using pap_solver::Algorithm;

typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS> Graph;

template<typename Graph>
static void parse_matrix(Graph* graph, std::istream* is) {
  static MatrixParser parser;
  auto add_edge_function = [](const typename Graph::vertex_descriptor& v1,
                              const typename Graph::vertex_descriptor& v2,
                              Graph* g) {
    boost::add_edge(v1, v2, *g);
  };
  parser.parse_full_matrix(is, graph, add_edge_function);
}

int main(int argc, char *argv[]) {
  Graph graph;
  std::ifstream file;
  file.open(argv[1]);
  parse_matrix(&graph, &file);

  Algorithm<Graph, std::default_random_engine> alg;
  alg.debug(graph);
  return 0;
}
