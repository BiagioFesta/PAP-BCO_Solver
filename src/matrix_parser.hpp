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


#ifndef __PAP_BCO_SOLVER__MATRIX_PARSER__HPP
#define __PAP_BCO_SOLVER__MATRIX_PARSER__HPP

#include <istream>
#include <ostream>
#include <vector>
#include <boost/graph/random.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>

namespace pap_solver {

class MatrixParser {
 public:
  template<typename Graph, typename RndEngine>
  static void generate_random_graph(const size_t num_vertices,
                                    const size_t num_edges,
                                    Graph* output_graph,
                                    RndEngine* rnd_engine);
  template<typename Graph>
  static void print_graph_raw(const Graph& graph, std::ostream* out_stream);

  template<typename Graph>
  static void read_graph_raw(std::istream* in_stream, Graph* graph);
};

template<typename Graph, typename RndEngine>
void MatrixParser::generate_random_graph(const size_t num_vertices,
                                         const size_t num_edges,
                                         Graph* output_graph,
                                         RndEngine* rnd_engine) {
  assert(output_graph != nullptr);
  assert(rnd_engine != nullptr);

  static_assert(std::is_integral<typename Graph::vertex_descriptor>::value,
                "The vertex type (descriptor) must to be a integer type!");

  output_graph->clear();

  boost::generate_random_graph(*output_graph,
                               num_vertices,
                               num_edges,
                               *rnd_engine,
                               false);

  const auto vertices = boost::vertices(*output_graph);
  int added_edges = 0;
  std::for_each(vertices.first,
                vertices.second,
                [&] (const typename Graph::vertex_descriptor& v) {
                  if (boost::out_degree(v, *output_graph) == 0) {
                    std::uniform_int_distribution<size_t> ovr(0, num_vertices);

                    auto ov = ovr(*rnd_engine);
                    while (ov == v) ov = ovr(*rnd_engine);

                    boost::add_edge(v, ov, *output_graph);
                    ++added_edges;
                  }
                });
}

template<typename Graph>
void MatrixParser::print_graph_raw(const Graph& graph,
                                   std::ostream* out_stream) {
  assert(out_stream != nullptr);

  boost::archive::text_oarchive archive(*out_stream);
  archive << graph;
}

template<typename Graph>
void MatrixParser::read_graph_raw(std::istream* in_stream,
                                  Graph* graph) {
  assert(in_stream != nullptr);
  assert(graph != nullptr);

  boost::archive::text_iarchive archive(*in_stream);
  archive >> *graph;
}


}  // namespace pap_solver

#endif
