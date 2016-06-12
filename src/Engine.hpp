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

  void find_a_solution(const Graph& graph,
                       Solution* output_solution,
                       int seed = -1);

  void find_a_solution_and_print(const Graph& graph,
                                 std::ostream* os,
                                 int seed = -1);

  bool check_solution(const Graph& graph,
                      const Solution& solution) const;
};


template<typename Graph>
void Engine<Graph>::find_a_solution_and_print(const Graph& graph,
                                              std::ostream* os,
                                              int seed) {
  AlgorithmDefault algorithm;
  if (seed >= 0) {
    algorithm.set_seed(seed);
  } else {
    algorithm.set_seed(
        std::chrono::system_clock::now().time_since_epoch().count());
  }

  Solution solution;
  algorithm.find_rnd_solution_fast(graph, &solution);

  *os << "--------Spanning Tree Considered-----------\n";
  solution.m_spanning_tree.print(graph, os);
  *os << "-------------------------------------------\n";
  *os << "----------------Solution-------------------\n";
  AlgorithmDefault::print_solution(os, solution);
  *os << "Number of vertices on PortAB: " <<
      solution.m_size_solution << "\n";
  *os << "Number of odd co-tree edges: " <<
      solution.m_number_odd_edges << "\n";
  *os << "Time elapsed for the solution: " <<
      solution.m_time_for_solution.count() << " ms\n";
  *os << "--------------------------------------------\n";
  assert(check_solution(graph, solution) == true);
}

template<typename Graph>
bool Engine<Graph>::check_solution(const Graph& graph,
                                   const Solution& solution) const {
  typedef typename AlgorithmDefault::VertexType VertexType;
  typedef typename AlgorithmDefault::PortAssignment PortAssignment;

  auto vertx = vertices(graph);
  auto finder = std::find_if(
      vertx.first,
      vertx.second,
      [&solution, &graph]
      (const VertexType& v) {
        // Get the assignment of v
        const auto& v_assignment = solution.m_assignment.at(v);
        if (v_assignment == PortAssignment::UnAssigned) {
          std::runtime_error("The solution is not complete!");
        } else if (v_assignment == PortAssignment::PortAB) {
          // If the node is PortAB it cannot give problem, so skip it
          return false;
        }

        // Get all adjacent and check if exsists one with
        // same assignment.
        auto adj_list = adjacent_vertices(v, graph);
        auto finder = std::find_if(
            adj_list.first,
            adj_list.second,
            [&solution, &v_assignment]
            (const VertexType& v_adj) {
              const auto& adj_assignment =
              solution.m_assignment.at(v_adj);

              if (adj_assignment == PortAssignment::UnAssigned) {
                std::runtime_error("The solution is not complete!");
              }

              if (adj_assignment == v_assignment) {
                return true;
              }
              return false;
            });

        if (finder != adj_list.second) {
          return true;
        }
        return false;
      });
  if (finder != vertx.second) {
    return false;
  }
  return true;
}
template<typename Graph>
void Engine<Graph>::find_a_solution(const Graph& graph,
                                    Solution* output_solution,
                                    int seed) {
  AlgorithmDefault algorithm;
  if (seed >= 0) {
    algorithm.set_seed(seed);
  } else {
    algorithm.set_seed(
        std::chrono::system_clock::now().time_since_epoch().count());
  }

  Solution& solution = *output_solution;
  algorithm.solve(graph, &solution);
  assert(check_solution(graph, solution) == true);
}


}  // namespace pap_solver


#endif  // __PAP_BCO_PARSER__ENGINE__HPP

