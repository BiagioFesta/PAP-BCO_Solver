#!/bin/bash
./create_datafile.py "$@" > data.dat && \
echo OPTIMAL SOLUTION && \
(time ampl ampl.run) 2>&1 | grep -P '(?:objective|real)' && \
echo PAP-SOLVER && \
(time ../../build/pap-bco_solver matrix.dat) 2>&1 | grep -P '(?:PortAB|real|edges)' && \
rm matrix.dat data.dat
