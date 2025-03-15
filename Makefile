# Compiler
CXX = g++

# ROOT flags
ROOT_FLAGS = $(shell root-config --cflags --glibs)

# Source file
SRC = VisEFiel.cpp

# Output executable
TARGET = VisEFiel

# Compilation rule
$(TARGET): $(SRC)
	$(CXX) -o $(TARGET) $(SRC) $(ROOT_FLAGS)

# Clean rule
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: clean
