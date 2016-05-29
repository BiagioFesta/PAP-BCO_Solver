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

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include "pap-bco_solver.hpp"
#include "matrix_parser.hpp"

int main(int argc, char *argv[]) {
  try {
    pap_solver::PAP_BCO_Solver solver(argc, argv);
    solver.run();
  } catch(const std::exception& err) {
    std::cerr << err.what() << '\n';
    return -1;
  }
  return 0;
}




namespace pap_solver {

PAP_BCO_Solver::PAP_BCO_Solver(int argc, char* argv[]) {
  if (parse_cmdline_options(argc, argv)) {
    throw std::invalid_argument("Invalid process"
                                " parsing command line arguments");
  }
}

int PAP_BCO_Solver::parse_cmdline_options(int argc, char* argv[]) {
  static struct option long_options[] = {
    {"file", required_argument, 0, 'f'},
    {"compressed", no_argument, 0, 'c'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };
  int option_index;
  while (true) {
    option_index = 0;
    auto c = getopt_long(argc, argv, "hcf:", long_options, &option_index);
    if (c == -1) break;
    switch (c) {
      case 'f':
        m_options.input_matrix_filename = optarg;
        break;
      case 'c':
        m_options.compressed_matrix = true;
        break;
      case 'h':
        m_options.display_help = true;
        break;
      case '?':
        return -1;
      default:
        std::cerr << "Invalid parsing options\n";
        return -1;
    }
  }
  if (m_options.display_help == false &&
      m_options.input_matrix_filename.size() == 0) {
    std::cerr << "You must specify an input file with the option: \n"
        "      -f FILENAME \n";
    return -1;
  }
  return 0;
}

void PAP_BCO_Solver::run() {
  if (m_options.display_help == true) {
    print_help();
    return;
  }

  parse_matrix_fromfile();
  MatrixParser parser;
  parser.print_graph_humanreadable(&m_graph, &std::cout);
}

void PAP_BCO_Solver::parse_matrix_fromfile() {
  MatrixParser parser;
  std::ifstream file;
  file.open(m_options.input_matrix_filename);
  if (file.fail() == true) {
    throw std::invalid_argument("The file '" +
                                m_options.input_matrix_filename +
                                "' cannot be read.");
  }
  if (m_options.compressed_matrix == true) {
    parser.parse_compressed_matrix(&file, &m_graph, AddEdge());
  } else {
    try {
      parser.parse_full_matrix(&file, &m_graph, AddEdge());
    } catch(const std::exception& err) {
      // TODO(biagio): make the suggestion more usefull
      std::string messages = err.what();
      messages += "\n\tMaybe the option '-c' could be usefull."
          " Use '-h' for help";
      throw std::runtime_error(messages);
    }
  }
  file.close();
}

void PAP_BCO_Solver::print_help() const noexcept {
  std::cout <<
      R"##(
   PAP-BCO Solver
-f, --file FILENAME     : The input matrix to load.
-c, --compresed         : Specify whether the input matrix is a compressed format or not.

)##";
}

}  // namespace pap_solver
