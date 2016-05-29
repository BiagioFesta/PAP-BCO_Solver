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
#include <stdexcept>
#include <memory>
#include <boost/graph/adjacency_list.hpp>
#include "matrix_parser.hpp"

namespace pap_solver {

template<typename Graph,
         typename AddEdge>
void MatrixParser::parse_compressed_matrix(std::istream* is,
                                           Graph* g,
                                           AddEdge add_edge) {
  // TODO(biagio): you should check whether the graph is empty or not.
  std::string sstream((std::istreambuf_iterator<char>(*is)),
                      std::istreambuf_iterator<char>());
  size_t current_vertex = 0;
  size_t current_column = 1;
  size_t prev_lenght_line = 0;
  for (const auto&c : sstream) {
    if (prev_lenght_line > 0 && current_column > prev_lenght_line) {
      throw std::runtime_error("Parsing error, line too long.");
    }
    switch (c) {
      case 48:  // ASCII for 0
        ++current_column;
        break;
      case 49:  // ASCII for 1
        add_edge(current_vertex, current_column, g);
        ++current_column;
        break;
      case 10:  // ASCII for \n
        if (prev_lenght_line > 0 && current_column != prev_lenght_line) {
          throw std::runtime_error("Parsing error, line too short.");
        }
        prev_lenght_line = current_column;
        ++current_vertex;
        current_column = current_vertex + 1;
      default:  // Pass away for other any char
        break;
    }
  }
}

template<typename Graph,
         typename AddEdge>
void MatrixParser::parse_full_matrix(std::istream* is,
                                     Graph* g,
                                     AddEdge add_edge) {
  static constexpr size_t SIZE_TEMP_BUFFER = 1024;
  auto tbuffer = std::get_temporary_buffer<char>(SIZE_TEMP_BUFFER);
  if (tbuffer.second != SIZE_TEMP_BUFFER) {
    throw std::runtime_error("Bad memory allocation");
  }

  size_t current_vertex = 0;
  size_t current_column = 1;
  size_t matrix_size = 0;
  size_t panning = 0;

  while (is->eof() == false) {
    panning = current_vertex + 1;
    is->getline(tbuffer.first, tbuffer.second);

    char* i = tbuffer.first;
    while (*i != 0) {
      switch (*i) {
        case 48:  // ASCII for 0
          if (panning > 0) {
            --panning;
          } else {
            ++current_column;
          }
          break;
        case 49:  // ASCII for 1
          if (panning > 0) {
            --panning;
          } else {
            add_edge(current_vertex, current_column, g);
            ++current_column;
          }
          break;
        default:  // Pass away for other any char
          break;
      }
      ++i;
    }
    if (matrix_size == 0) {
      matrix_size = current_column;
    } else {
      if (matrix_size != current_column) {
        throw std::runtime_error("Parsing error");
      }
    }
    ++current_vertex;
    current_column = current_vertex + 1;
  }

  std::return_temporary_buffer(tbuffer.first);
}


template<typename Graph>
void MatrixParser::print_graph_humanreadable(Graph* g,
                                             std::ostream* os) {
  auto its = boost::vertices(*g);
  for (auto i = its.first; i != its.second; ++i) {
    auto adjacent = boost::adjacent_vertices(*i, *g);
    for (auto j = adjacent.first; j != adjacent.second; ++j) {
      *os << *i << "  ->  " << *j << '\n';
    }
  }
}

}  // namespace pap_solver
