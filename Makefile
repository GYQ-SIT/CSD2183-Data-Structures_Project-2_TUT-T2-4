CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

all: simplify validate

simplify: main.cpp apsc.h ring.h geometry.h spatial_index.h
	$(CXX) $(CXXFLAGS) -o simplify main.cpp

validate: validate.cpp
	$(CXX) $(CXXFLAGS) -o validate validate.cpp

clean:
	rm -f simplify validate
