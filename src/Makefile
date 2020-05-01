CXX = g++
CXXFLAGS = -std=c++0x -O3 -Wall -fopenmp -Wno-maybe-uninitialized
CXXFILES = main.cpp

all: a.out

a.out: $(CXXFILES)
		$(CXX) $(CXXFLAGS) -o a.out $(CXXFILES)

clean:
		rm a.out
