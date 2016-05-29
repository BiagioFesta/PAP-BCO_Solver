CXX = g++
CXXFLAGS = -g -Wall
OUTPUT = pap-bco_solver
SRC = src/pap-bco_solver.cpp

all:
	${CXX} ${CXXFLAGS} -o ${OUTPUT} ${SRC}
