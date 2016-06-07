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

#ifndef __PAP_BCO_PARSER__GRAPH_UTILITY__HPP
#define __PAP_BCO_PARSER__GRAPH_UTILITY__HPP

#include <queue>
#include <boost/graph/two_bit_color_map.hpp>

namespace pap_solver {

class GraphUtility {
 public:
  template<typename Graph>
  static bool graph_contains_loop(const Graph& graph);
};

template<typename Graph>
bool GraphUtility::graph_contains_loop(const Graph& graph) {
  // Type definitions
  typedef typename Graph::vertex_descriptor VertexType;
  typedef typename Graph::edge_descriptor EdgeType;
  typedef boost::two_bit_color_map<> color_map_t;
  typedef std::queue<VertexType> open_list_t;

  // Color definitions
  static const auto white_t =  boost::color_traits<
    boost::two_bit_color_type>::white();
  static const auto black_t =  boost::color_traits<
    boost::two_bit_color_type>::black();
  static const auto gray_t =  boost::color_traits<
    boost::two_bit_color_type>::gray();

  const auto num_vertices = boost::num_vertices(graph);
  color_map_t color_map(num_vertices);

  open_list_t openlist;
  const auto root = *(vertices(graph).first);
  openlist.push(root);
  put(color_map, root, black_t);

  while (openlist.empty() == false) {
    const auto& node = openlist.front();
    auto this_node_color = get(color_map, node);
    auto edges_list = out_edges(node, graph);
    for (auto i = edges_list.first;
         i != edges_list.second;
         ++i) {
      auto target = boost::target(*i, graph);
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

}  // namespace pap_solver

#endif  // __PAP_BCO_PARSER__GRAPH_UTILITY__HPP
