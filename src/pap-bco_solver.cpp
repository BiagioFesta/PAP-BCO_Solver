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
#include <string>
#include <fstream>
#include <stdexcept>
#include <ctime>
#include <queue>
#include <map>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include "pap-bco_solver.hpp"
#include "Engine.hpp"

int main(int argc, char *argv[]) {
  try {
    pap_solver::PAP_BCO_Solver solver;
    solver.run(argc, argv);
  } catch(const std::exception& err) {
    std::cerr << err.what() << '\n';
    return -1;
  }
  return 0;
}




namespace pap_solver {

PAP_BCO_Solver::PAP_BCO_Solver() noexcept {
}

int PAP_BCO_Solver::parse_cmdline_options(int argc, char* argv[]) {
  static struct option long_options[] = {
    {"compressed", no_argument, 0, 'c'},
    {"help", no_argument, 0, 'h'},
    {"generate", required_argument, 0, 'g'},
    {"seed", required_argument, 0, 's'},
    {0, 0, 0, 0}
  };
  int option_index;
  while (true) {
    option_index = 0;
    auto c = getopt_long(argc, argv, "s:g:hc", long_options, &option_index);
    if (c == -1) break;
    switch (c) {
      case 'c':
        m_options.compressed_matrix = true;
        break;
      case 'h':
        m_options.display_help = true;
        break;
      case 's':
        m_options.debug_seed = true;
        try {
          m_options.seed = std::stoi(optarg);
        } catch(const std::invalid_argument& err) {
          std::cerr << "--seed [NUMBER]   Number must to be an"
              " integer number";
          return -1;
        }
        break;
      case 'g': try {
          std::string temp_arg = optarg;
          auto finder = temp_arg.find(':');
          if (finder == std::string::npos) {
            throw std::invalid_argument("NEDGES not specified.");
          }

          m_options.generate_random_graph = true;

          m_options.generate_num_vertices =
              std::stoi(temp_arg.substr(0, finder));

          m_options.generate_num_edges =
                std::stoi(temp_arg.substr(finder + 1));
        } catch (const std::exception& err) {
          std::cerr << "--generate NVERTICES:NEDGES\n";
          return -1;
        }
        break;
      case '?':
        return -1;
      default:
        std::cerr << "Invalid parsing options\n";
        return -1;
    }
  }
  if (optind < argc) {
    m_options.input_filename = argv[optind];
  }
  return 0;
}

void PAP_BCO_Solver::run(int argc, char* argv[]) {
  if (parse_cmdline_options(argc, argv)) {
    throw std::invalid_argument(
        "Invalid process parsing command line arguments");
  }

  if (m_options.display_help == true) {
    print_help();
  } else {
    if (m_options.generate_random_graph == true) {
      auto seed = m_options.debug_seed ? m_options.seed : std::time(nullptr);
      std::default_random_engine rnd_eng(seed);

      generate_random_matrix(&rnd_eng);
    } else {
      print_header();

      parse_matrix();
      std::cout << "------------Matrix Information-------------\n";
      std::cout << "Number of vertices: " << boost::num_vertices(m_graph) <<
          "\nNumber of edges: " << boost::num_edges(m_graph) << "\n";
      std::cout << "-------------------------------------------\n";

      Engine<Graph> engine_solver;
      engine_solver.find_a_solution_and_print(
          m_graph,
          &std::cout,
          m_options.debug_seed ? m_options.seed : -1);
    }
  }
}

void PAP_BCO_Solver::parse_matrix() {
  std::istream* is = &std::cin;
  std::ifstream file;
  if (m_options.input_filename.size() != 0) {
    file.open(m_options.input_filename);
    if (file.fail() == true)
      throw std::invalid_argument(
          "The file '" + m_options.input_filename + "' cannot be open.");
    is = &file;
  }

  GraphUtility::readGraph_asArchive(is, &m_graph);

  if (file.is_open()) {
    file.close();
  }
}

void PAP_BCO_Solver::print_help() const noexcept {
  std::cout <<
      R"##(

Use:    pap-bco_solver [OPTION]... [FILENAME]
   Command Options:

FILENAME                       It's the input (or output if --generate option has been
                               specified) of the adjacency matrix.
                               If it's not specified then the matrix will be read from
                               the std input (or generated to the std output).

-c, --compresed                Specify whether the input matrix is a compressed format or not.

-g NVERTICES:NEDGES,
--generate NVERTICES:NEDGES    Generate a valid random graph  with NVERTICES and NEDGES.
                               Both those numbers must to be integer numbers.

-h, --help                     Display this guide.

-s, --seed NUM                 Specify the seed for random generator.

)##";
}

void PAP_BCO_Solver::print_header() const noexcept {
  std::cout << R"##(
           PAP-BCO Solver
   Copyright (C) 2016  Biagio Festa
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

)##";
}

}  // namespace pap_solver
