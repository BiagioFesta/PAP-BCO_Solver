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

#ifndef __PAP_BCO_PARSER__ENGINE__HPP
#define __PAP_BCO_PARSER__ENGINE__HPP

#include <random>
#include "Algorithm.hpp"

namespace pap_solver {

template<typename Graph>
class Engine {
 public:
  typedef std::default_random_engine RndGenerator;
  typedef Algorithm<Graph, RndGenerator> AlgorithmDefault;
  typedef typename AlgorithmDefault::Solution Solution;

  Engine() = default;
  void find_a_solution_and_print(const Graph& graph,
                                 std::ostream* os);
};


template<typename Graph>
void Engine<Graph>::find_a_solution_and_print(const Graph& graph,
                                              std::ostream* os) {
  AlgorithmDefault algorithm;
#ifdef _DEBUG
  algorithm.set_seed(0);
#else
  algorithm.set_seed(
      std::chrono::system_clock::now().time_since_epoch().count());
#endif

  Solution solution;
  algorithm.solve(graph, &solution);

  *os << "--------Spanning Tree Considered-----------\n";
  // TODO(biagio): print spanning tree
  *os << "-------------------------------------------\n";
  *os << "----------------Solution-------------------\n";
  AlgorithmDefault::print_solution(os, solution);
  *os << "Number of vertices on PortAB: " <<
      solution.m_size_solution << "\n";
  *os << "Time elapsed for the solution: " <<
      solution.m_time_for_solution.count() << " ms\n";
  *os << "--------------------------------------------\n";
}

}  // namespace pap_solver


#endif  // __PAP_BCO_PARSER__ENGINE__HPP

