# Compiler
CXX = g++

# Garfield++ and ROOT flags
GARFIELD_INCLUDE = -I$(GARFIELD_HOME)/install/include/
GARFIELD_LIBS = -L$(GARFIELD_HOME)/install/lib/ -lGarfield
ROOT_FLAGS = $(shell root-config --cflags --glibs)

# Source file
SRC = GALASim.cpp

# Output executable
TARGET = GALASim

# Compilation rule
$(TARGET): $(SRC)
	$(CXX) -o $(TARGET) $(SRC) $(GARFIELD_INCLUDE) $(GARFIELD_LIBS) $(ROOT_FLAGS)

# Clean rule
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: clean
