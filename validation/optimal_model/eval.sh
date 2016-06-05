#!/bin/bash

./create_datafile_from_matrix.py "$@" > data.dat && \
echo OPTIMAL SOLUTION && \
(time ampl ampl.run) 2>&1 | grep -P '(?:objective|real)' && \
echo PAP-SOLVER && \
(time ../../build/pap-bco_solver $1) 2>&1 | grep -P '(?:PortAB|real|edges)' && \
rm data.dat
