all: compile execute

compile:
	g++ -std=c++11 -O3 -flto -Wall -o CodeBase/main.x CodeBase/main.cpp CodeBase/FluidDynamics.cpp  -larmadillo -llapack -lblas
execute:
	./main.x
