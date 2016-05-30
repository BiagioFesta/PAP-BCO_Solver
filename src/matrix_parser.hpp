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
  void parse_compressed_matrix(std::istream* is, Graph* g, AddEdge add_edge);

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
  void parse_full_matrix(std::istream* is, Graph* g, AddEdge add_edge);


  /// @brief Generate a random compressed matrix
  ///
  /// @param [in] rnd_engine     A valid random engine.
  /// @param [in] matrix_size    The size of the matrix.
  /// @param [out] os            An output stream where the matrix
  ///                            will be written.
  template<typename RND>
  void generate_rnd_compressed_matrix(RND* rnd_engine,
                                      size_t matrix_size,
                                      std::ostream* os);

  template<typename Graph>
  void print_graph_humanreadable(Graph* g, std::ostream* os);
};

}  // namespace pap_solver

#include "matrix_parser.t.hpp"
#endif
