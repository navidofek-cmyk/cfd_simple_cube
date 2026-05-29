CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++20 -Wall -Iinclude -fopenmp
SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:.cpp=.o)
BIN := channelCube

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(BIN)
	./$(BIN)

clean:
	rm -f $(OBJ) $(BIN) *.dat *.vtk

.PHONY: run clean
