CXX ?= g++

BOOST_INCLUDE ?= /tmp/boost_1_83_0
WSTP_DIR ?= /usr/local/Wolfram/Mathematica/14.3/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions

CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp
CXXFLAGS += -I$(BOOST_INCLUDE) -I$(WSTP_DIR)
LDLIBS := -L$(WSTP_DIR) -Wl,-rpath,$(WSTP_DIR) -lWSTP64i4 -lmpc -lmpfr -lgmp -lpthread -lrt -luuid -ldl -fopenmp

TARGET := chess_native_link
GENERATED := src/linktm.cpp
OBJECTS := build/backend.o build/linktm.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

$(GENERATED): src/link.tm
	$(WSTP_DIR)/wsprep $< -o $@

build/backend.o: src/backend.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/linktm.o: $(GENERATED)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build $(TARGET) $(GENERATED)
