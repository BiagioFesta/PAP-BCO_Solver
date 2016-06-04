FROM ubuntu
MAINTAINER biagio.festa[at]gmail.com

RUN apt-get update && \
	apt-get install -y gcc g++ git cmake libboost-all-dev

WORKDIR ~/

RUN git clone https://github.com/BiagioFesta/PAP-BCO_Solver.git && \
	cd PAP-BCO_Solver && \
	cmake . && \
	make                                                                    