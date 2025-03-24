# Compiler
CXX = g++

# Garfield++ and ROOT flags
GARFIELD_INCLUDE = -I$(GARFIELD_HOME)/install/include/
GARFIELD_LIBS = -L$(GARFIELD_HOME)/install/lib/ -lGarfield
ROOT_FLAGS = $(shell root-config --cflags --glibs)

# Source files and output executables
SRC1 = GALASim.cpp
SRC2 = GALASim_clu.cpp
TARGET1 = GalaSim
TARGET2 = GALASim_clu

# Compilation rules
all: $(TARGET1) $(TARGET2)

$(TARGET1): $(SRC1)
	$(CXX) -o $(TARGET1) $(SRC1) $(GARFIELD_INCLUDE) $(GARFIELD_LIBS) $(ROOT_FLAGS)

$(TARGET2): $(SRC2)
	$(CXX) -o $(TARGET2) $(SRC2) $(GARFIELD_INCLUDE) $(GARFIELD_LIBS) $(ROOT_FLAGS)

# Clean rule
clean:
	rm -f $(TARGET1) $(TARGET2)

# Phony targets
.PHONY: all clean
