// Copyright 2016 <Biagio Festa>

#include <string>
#include <memory>
#include <cstring>
#include "../src/matrix_parser.hpp"
#include "../src/Engine.hpp"
using pap_solver::MatrixParser;
using pap_solver::Engine;
using pap_solver::SpanningTree;

typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS> Graph;

static const std::string embedded_matrix =
    "110\n"
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
  decltype(engine)::Solution solution;
  parse_matrix(&graph);
  engine.find_a_solution(graph, &solution, 0);

  const auto e02 = edge(0, 2, graph);
  const auto e13 = edge(1, 3, graph);
  assert(e02.second == true);
  assert(e13.second == true);
  
  solution.m_spanning_tree.print(graph, &std::cout);
  solution.m_spanning_tree.perform_transformation(graph,
                                                  e13.first,
                                                  e02.first);
  

  return 0;
}
