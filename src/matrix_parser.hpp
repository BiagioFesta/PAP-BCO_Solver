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
#include <string>

namespace pap_solver {

class MatrixParser {
 public:
  /// \brief It parses a compressed matrix and fills the relative graph.
  ///        The matrix must to be in the compressed form:
  ///              a12   a13    a14    ...   a1n
  ///              a23   a24    ...    a2n
  ///              ...
  ///              a(n-1)(n)
  ///
  /// \param [in] is        An input stream where the matrix is stored
  ///                       as text file.
  /// \param [out] g        The graph has to be generated.
  ///                       The graph passed must to be empty.
  /// \param [in] add_edge  Is a object function which performs the operation
  ///                       of insertion of an edge into the graph.
  ///                       Concept:    void (&)(VertexType, VertexType, Graph)
  template<typename Graph,
           typename AddEdge>
  void parse_compressed_matrix(std::istream* is,
                               Graph* g,
                               AddEdge add_edge) const;

  /// \brief It prarses a full matrix (not compressed) and filles the relative
  ///        graph
  /// \param [in] is       An input stream where the matric is stored
  ///                      as text file.
  /// \param [out] g       The graph has to be generated.
  ///                      The graph passed must to be empty.
  /// \param [in] add_edge Is an object function which perform the operatiorn
  ///                      of insertion of an edge into the graph.
  /// \see                 MatrxiParser::parse_compressed_matrix.
  template<typename Graph,
           typename AddEdge>
  void parse_full_matrix(std::istream* is,
                         Graph* g,
                         AddEdge add_edge) const;


  /// @brief Generate a random compressed matrix
  ///
  /// @param [in] rnd_engine     A valid random engine.
  /// @param [in] matrix_size    The size of the matrix.
  /// @param [out] os            An output stream where the matrix
  ///                            will be written.
  template<typename RND>
  void generate_rnd_compressed_matrix(RND* rnd_engine,
                                      size_t matrix_size,
                                      std::ostream* os) const;

  template<typename RND>
  void generate_rnd_matrix(RND* rnd_engine,
                           size_t matrix_size,
                           std::ostream* os) const;
};




template<typename Graph,
         typename AddEdge>
void MatrixParser::parse_compressed_matrix(std::istream* is,
                                           Graph* g,
                                           AddEdge add_edge) const {
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
        } else if (prev_lenght_line == 0 && current_column == 1) {
          throw std::runtime_error("Cannot parse an empty matrix");
        }
        prev_lenght_line = current_column;
        ++current_vertex;
        current_column = current_vertex + 1;
      default:  // Pass away for other any char
        break;
    }
  }
  if (current_vertex == 0) {
    throw std::runtime_error("Cannot parse an empty matrix");
  }
}

template<typename Graph,
         typename AddEdge>
void MatrixParser::parse_full_matrix(std::istream* is,
                                     Graph* g,
                                     AddEdge add_edge) const {
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
    if (std::strlen(tbuffer.first) == 0) continue;
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
  if (current_vertex == 0) {
    throw std::runtime_error("Cannot parse an empty matrix");
  }
  std::return_temporary_buffer(tbuffer.first);
}

template<typename RND>
void MatrixParser::generate_rnd_compressed_matrix(RND* rnd_engine,
                                                  size_t matrix_size,
                                                  std::ostream* os) const {
  // 48 and 49 are '0' and '1' in the ASCII codec
  std::uniform_int_distribution<> rnd_value(48, 49);

  auto tbuffer = std::get_temporary_buffer<char>(matrix_size);
  if (tbuffer.first == nullptr) {
    throw std::runtime_error("I cannot generate random matrix, "
                             "bad memory allocation");
  }
  size_t i;
  size_t row_size = matrix_size - 1;
  while (row_size) {
    for (i = 0; i < row_size; ++i) {
      tbuffer.first[i] = rnd_value(*rnd_engine);
    }
    // We cannot have a vertex without outgoing edge
    if (std::memchr(tbuffer.first, '1', row_size) == nullptr) {
      // We have a row with all zeros. Just change on of them!
      tbuffer.first[0] = '1';
    }
    os->write(tbuffer.first, row_size);
    os->put('\n');
    --row_size;
  }
  std::return_temporary_buffer(tbuffer.first);
}

template<typename RND>
void MatrixParser::generate_rnd_matrix(RND* rnd_engine,
                                       size_t matrix_size,
                                       std::ostream* os) const {
  if (matrix_size == 0) return;
  std::uniform_int_distribution<> rnd_value(48, 49);
  char* matrix_data = new char[matrix_size*matrix_size];
  for (auto i = 0; i < matrix_size; ++i) {
    // Generate the i-th row of the matrix
    if (i > 0) {
    }
    for (auto j = 0; j < i; ++j) {
      // Generate the j-th column of the matrix with known data.
      matrix_data[i*matrix_size + j] = matrix_data[j*matrix_size + i];
    }
    for (auto j = i; j < matrix_size; ++j) {
      // Generate the j-th column of the matrix with new data.
      char data_element;
      if (j == i) {
        data_element = '0';
      } else {
        data_element = static_cast<char>(rnd_value(*rnd_engine));
      }
      matrix_data[i*matrix_size + j] = data_element;
    }
    if (std::memchr(matrix_data + (i*matrix_size),
                    '1',
                    matrix_size) == nullptr) {
      // We have a row with all zeros. Just change on of them!
      matrix_data[i*matrix_size + matrix_size - 1] = '1';
    }
    os->write(matrix_data + (i*matrix_size), matrix_size);
    os->put('\n');
  }
  delete[] matrix_data;
}

}  // namespace pap_solver

#endif
