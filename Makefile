INCFLAGS = -I/usr/local/include/ -I./util/

CPP = mpicxx
CPPFLAGS = -g -O3 $(INCFLAGS)  -fopenmp -Wall -Wno-strict-aliasing -std=c++11
LINKERFLAGS = -lz -ltbb
DEBUGFLAGS = -g -ggdb $(INCFLAGS)
HEADERS=$(shell find . -name '*.hpp')


all: src/GraphSim_master src/GraphSim_worker src/GraphSim_swap_worker
echo:
	echo $(HEADERS)
clean:
	@rm -rf bin/*

src/% : src/%.cpp $(HEADERS)
	@mkdir -p bin/$(@D)
	$(CPP) $(CPPFLAGS) -Isrc/ $@.cpp -o bin/$@ $(LINKERFLAGS) 



