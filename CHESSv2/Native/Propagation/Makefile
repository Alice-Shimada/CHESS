CXX ?= g++

BOOST_INCLUDE ?= /tmp/boost_1_83_0
WSTP_DIR ?= /usr/local/Wolfram/Mathematica/14.3/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions
FLINT_DIR ?= /home/liuyuanche/MMA_Package/flint-3.4.0
FLINT_INCLUDE ?= /usr/local/include

CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp
CXXFLAGS += -I$(BOOST_INCLUDE) -I$(WSTP_DIR) -I$(FLINT_INCLUDE)
LDLIBS := -L$(WSTP_DIR) -Wl,-rpath,$(WSTP_DIR) -lWSTP64i4 -L$(FLINT_DIR) -Wl,-rpath,$(FLINT_DIR) -lflint -lmpc -lmpfr -lgmp -lpthread -lrt -luuid -ldl -fopenmp

TARGET := chess_native_link
GENERATED := src/linktm.cpp
OBJECTS := build/backend.o build/flint_coupled_solver.o build/linktm.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

$(GENERATED): src/link.tm
	$(WSTP_DIR)/wsprep $< -o $@

build/backend.o: src/backend.cpp src/backend_types.hpp src/flint_coupled_solver.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/flint_coupled_solver.o: src/flint_coupled_solver.cpp src/flint_coupled_solver.hpp src/backend_types.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/linktm.o: $(GENERATED)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build $(TARGET) $(GENERATED)
