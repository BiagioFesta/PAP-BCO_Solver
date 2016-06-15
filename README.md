# PAP-BCO Solver
Solver for "Port Assignment Problem for Binary Commutative Operators" (PAP-BCO)

## Introduction ##
Given a binding of operations and variables, the PAP-BCO consists in 
minimizing the number of interconnection between and registers and function unit
which implements that operation.

For more details about the problem and the proposed solution
see the slide-presentation at:

[http://www.biagiofesta.it/pap-bco_solver.php](http://www.biagiofesta.it/pap-bco_solver.php).


## Target Platform ##
The current platform supported is Linux.

## Installation ##
Use CMake:
* cd pap-bco_solver/
* cmake . -DCMAKE_BUILD_TYPE=Release
* make
* make install (optional)

