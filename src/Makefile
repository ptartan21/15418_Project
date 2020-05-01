CXX = g++
CXXFLAGS = -std=c++0x -Wall -fopenmp -Wno-maybe-uninitialized
CXXFILES = main_scc.cpp

all: a.out

a.out: $(CXXFILES)
		$(CXX) $(CXXFLAGS) -o a.out $(CXXFILES)

clean:
		rm a.out
