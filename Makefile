CC=mpicc
CXX=mpic++
ifeq ($(DEBUG),yes)
	CXXFLAGS=-std=c++11 -g -fopenmp -pthread -Wall -pedantic
else
	CXXFLAGS=-std=c++11 -O3 -fopenmp -pthread -Wall -pedantic
endif

CFLAGS=-std=c99 -O3 -fopenmp -pthread -Wall -Wextra

EXEC= raytracer.exe
SRC= raytracer.cpp
OBJ= $(SRC:.cpp=.cpp.o)

BIN = raytracer
CSRC = raytracer.c

$(BIN): $(CSRC)
	@$(CC) -o $@ $^ $(CFLAGS)


$(EXEC): $(SRC)
	@$(CXX) -o $@ $^ $(CXXFLAGS)


.PHONY: clean cleanall

clean:
	@rm -fr *.o *~

cleanall: clean
	@rm -fr $(EXEC)
	
