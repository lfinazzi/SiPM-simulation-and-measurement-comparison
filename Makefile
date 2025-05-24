# Compiler
CXX := g++

# Compiler flags
CXXFLAGS := -Wall -Wextra -std=c++17 $(shell root-config --cflags)

# Linker flags
LDFLAGS := $(shell root-config --libs)

# Target executable
TARGET := main

# Source files
SRCS := main.cpp aux_functions.cpp realdata.cpp simulation.cpp minimizer.cpp

# Object files
OBJS := $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean