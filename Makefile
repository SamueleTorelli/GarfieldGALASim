# Compiler
CXX = g++

# ROOT flags
ROOT_FLAGS = $(shell root-config --cflags --glibs)

# Source file
SRC = VisEFiel.cpp
SRC2 = EFSim3x3.cpp
SRC3 = EFSim3x3_SOR.cpp
SRC4 = EFSim1x1_SOR.cpp
SRC5 = EFSim1x1_SOR_2D.cpp

# Output executable
TARGET = VisEFiel
TARGET2 = EFSim3x3
TARGET3 = EFSim3x3_SOR
TARGET4 = EFSim1x1_SOR
TARGET5 = EFSim1x1_SOR_2D

# Compilation rules                                                                                                                                                                                                                                      
all: $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5)

# Compilation rule
$(TARGET): $(SRC)
	$(CXX) -o $(TARGET) $(SRC) $(ROOT_FLAGS)

$(TARGET2): $(SRC2)
	$(CXX) -o $(TARGET2) $(SRC2) $(ROOT_FLAGS)

$(TARGET3): $(SRC3)
	$(CXX) -o $(TARGET3) $(SRC3) $(ROOT_FLAGS)

$(TARGET4): $(SRC4)
	$(CXX) -o $(TARGET4) $(SRC4) $(ROOT_FLAGS)

$(TARGET5): $(SRC5)
	$(CXX) -o $(TARGET5) $(SRC5) $(ROOT_FLAGS)

# Clean rule
clean:
	rm -f $(TARGET) $(TARGET2) $(TARGET3) $(TARGET4)

# Phony targets
.PHONY: clean
